set -e

welcome(){
    help="
| ===========================================================================\n|\n\
| \t    Welcome to showrunner \n\
| \t'Together we run the show' ~ Kat DeLuna :) \n|\n\
| =========================================================================== \n|\n\
| Software used in this script requires Python > 3.8 \n\
| Please ensure they are installed \n\
| Recommend running this in a virtual environment \n|\n\
| --------------------------------------------------------------------------- \n\
| \tInstructions \n\
| \t1) Running everything at once  \n\
| ---------------------------------------------------------------------------\n\
\t1. Place this script in a sub-dir where your channelised images are \n\
  \ti.e. the channelised images should be at '../' with respect to where\n\
     \tthis script is \n\
\t2. Make script executable as: 'chmod +x showrunner.sh' \n\
\t3. Run the script as follows: \n\
\t\t ./run_most_rm_things location/of/mask/dir \n\
\te.g \n \
\t\t ./run_most_rm_things.sh /home/andati/pica/reduction/experiments/emancipation/masks-572\n\n\
| ---------------------------------------------------------------------------\n\
"

    echo -e $help

}

installRequiredSoftware(){
    # pip install -U pip
    pkgs=("MontagePy" "Owlcat" "spimple")

    echo -e "\n############################################################"
    echo "Installing required packages if need be"
    echo -e "\n############################################################"

    for pkg in ${pkgs[@]}
        do 
            if ! pip list --local | grep -i "$pkg";
            then
                echo "Package $pkg was not found, installing it."
                pip install $pkg
            else
                echo "$pkg already installed. Skipping it."
            fi
        done
    return 0
}


initialiseEnvVarsForThisScript(){
    echo -e "\n############################################################"
    echo "We initialise some of the environment variables in env-vars"
    echo -e "\n############################################################"
    if [[ ! -f env-vars ]]
    then
    	echo "env-vars does not exist. Creating"
    	ln -s $HOME/git_repos/misc_scripts_n_tools/env-vars
    fi

    source env-vars

    if [[ -z $mask_dir ]]
    then
        export mask_dir=masks

        if [[ ! -d $mask_dir ]]
        then
            mkdir -p $mask_dir
            echo -e "Created mask directory: $mask_dir"
        else
            echo -e "Mask directory already exists:  $mask_dir"
        fi
    fi

    return 0
}



makeDirs(){
    echo -e "\n############################################################"
    echo "Make these directories"
    echo -e "############################################################\n"
    mkdir -p $ORIG_CUBES $SEL_CUBES $CONV_CUBES $PRODS $SPIS $IMGS
    return 0
}



selectGoodChannels(){
    
    # 1. Create my channel selections here
    if [[ ! -f selected-channels.txt ]]
    then
    	echo -e "\n############################################################"
        echo -e """Autoselect some valid channels. Should be stored in a file 
                \r called selected-channels"""
        echo -e "############################################################\n"
        sc-beam-plot --output ./ --prefix 00 --threshold 0.5 --auto-select ../
    fi

    mapfile -t sel < selected-channels.txt; export sel;

    
    # 2. copy those selected channels' images to some other directory
    echo -e "\n############################################################"
    echo "copy relevant channels' images to this folder for IQUV"
    echo -e "############################################################\n"
    for n in ${sel[@]}
    	do
    		cp ../*-$n-{I,Q,U}-image.fits $IMGS
    	done

    
    # 3. generate frequency file for the selected images if its not there already
    if [[ ! -f frequencies.txt ]]
    then
        echo -e "\n############################################################"
        echo -e """write out the selected freqs into a single file for easier work. 
            \r This is for use in the RM synth for wavelength"""
        echo -e "############################################################\n"
        for im in $IMGS/*-I-image.fits
            do
                # ensure that this is being appended Lexy !!!!!!
                fitsheader -k CRVAL3 $im |  grep -i CRVAL3 >> frequencies.txt
            done


        echo -e "\n############################################################"
        echo -e """Cleanup freqs file by replacing CRVAL3 = and all spaces after it 
                \r with emptiness"""
        echo -e "############################################################\n"
        sed -i "s/CRVAL3  =\ \+//g" frequencies.txt
    fi


    
    # 4. save the names of the selected images in a file
    echo -e "\n############################################################"
    echo -e """Save the names of the selected images"""
    echo -e "############################################################\n"
    ls $IMGS/*-[0-9][0-9][0-9][0-9]*-image.fits > selected-freq-images.txt

    #Replacing all begins of strings here with ../
    sed -i 's/^/\.\.\//g' selected-freq-images.txt

    return 0
}


stackAllImages(){
    # Create an image cube from channelised images for Stokes IQU

    echo -e "\n############################################################"
    echo -e """Generating stokes cubes from channelised images"""
    echo -e "############################################################\n"

    for s in $STOKES
    do
        echo "make cubes from ALL the output images: ${s^^}"
        echo -e "---------------------------------------------------------------\n"
        
        images=$(ls -v ../*-[0-9][0-9][0-9][0-9]-$s-image.fits)
        fitstool.py --stack=$ORIG_CUBES/${s,,}-cube.fits:FREQ $(echo $images)

        if [[ ${s,,} = "i" ]]
        then
            images=$(ls -v ../*-[0-9][0-9][0-9][0-9]-$s-model.fits)
            fitstool.py --stack=$ORIG_CUBES/${s,,}-model-cube.fits:FREQ $(echo $images)
            images=$(ls -v ../*-[0-9][0-9][0-9][0-9]-$s-residual.fits)
            fitstool.py --stack=$ORIG_CUBES/${s,,}-residual-cube.fits:FREQ $(echo $images)
        fi
    done

    return 0
}


selectedChannels_Stack(){
    # optional argument syntax: variable=${read_arg_in_position:-default_value}
    local indir=${1:-$IMGS}
    local outdir=${2:-$SEL_CUBES}
    local suffix=${3:-"image-cube"}

    # print function information
    echo -e "############################################################\n"
    echo "Running function: $FUNCNAME($@)"
    echo -e "############################################################\n"
    
    for s in $STOKES
        do
            echo "Make the selection cubes: ${s^^}"
            images=$(ls -v $indir/*-$s-image.fits)
            fitstool.py --stack=$outdir/${s,,}-$suffix.fits:FREQ $(echo $images)
            echo "Stored at: $outdir"
        done
    return 0
}


selectedChannels_CubeConvolve(){    
    # 1. Stack all the SELECTED channel images to form image cubes
    # 2. Convolve the CUBES to the same resolution

    mapfile -t beam_dims < beam-dims.txt
    
    echo -e "\n############################################################"
    echo "Convolve the cubes to the same resolution"
    echo -e "############################################################\n"
    for s in $STOKES
        do        
            spimple-imconv -image $SEL_CUBES/${s,,}-image-cube.fits \
                -o $CONV_CUBES/${s,,}-image-cube -pp ${beam_dims[@]}
        done


    echo -e "\n############################################################"
    echo "Renaming output file from spimple because the naming here is weird"
    echo -e "############################################################\n"
    rename.ul -- ".convolved.fits" ".fits" $CONV_CUBES/* || true


    echo -e "\n############################################################"
    echo "Just check if the beam sizes are the same"
    echo -e "############################################################\n"
    # sc-beam-plot --output ./ --threshold 0.5  $IMGS

    return 0

}



selectedChannels_SinglesConvolve(){

    # 1. Convolve each selected channel and store it in a different directory
    # 2. Stack all the convolved channel images to form image cubes
    # 3. Plot the beam sizes of the new convolved channelised images

    mkdir -p $CONVIM

    echo -e "\n############################################################"
    echo -e """Convolve each channelised image INDIVIDUALLY"""
    echo -e "############################################################\n"

    if [[ ! -f beam-dims.txt ]]
    then
        # read beam dimensions of the first available channel
        fitsheader $IMGS/$(ls -v relevant-images/ | head -1) \
        |   grep -i "bmaj \|bmin \|bpa " > beam-dims.txt
        sed -i "s/.*=\s*//g" beam-dims.txt

    fi
    # beam_dims=$(python -c "import numpy as np; cat = np.loadtxt('beam-dims.txt'); print(*cat)")
    mapfile -t beam_dims < beam-dims.txt

    for im in $(ls -v $IMGS/*.fits)
    do
        spimple-imconv -image $im -pp ${beam_dims[@]} -o $CONVIM/$(basename $im)
    done

    rm $CONVIM/*.clean_psf*
    rename.ul ".convolved.fits" "" $CONVIM/*

    # stack them
    selectedChannels_Stack $CONVIM $CONV_CUBES conv-image-cube


    echo -e "\n############################################################"
    echo "Just check if the beam sizes are the same"
    echo -e "############################################################\n"
    sc-beam-plot --output ./ --prefix 01 $IMGS

    return 0
}


copyMfsImages(){
    echo -e "\n############################################################"
    echo "copy I MFS image reference image here"
    echo -e "############################################################\n"
    
    cp ../*MFS-I-image.fits i-mfs.fits
    cp ../*MFS-Q-image.fits q-mfs.fits
    cp ../*MFS-U-image.fits u-mfs.fits

}


generateSpiMap(){

    echo -e "\n############################################################"
    echo "Generating spectral index maps"
    echo -e "############################################################\n"
    
    # Generate some spectral index maps using spimple

    if [[ ! -f wsums.txt ]]
    then
        echo -e "\n############################################################"
        echo "Get their wsums and store"
        echo -e "############################################################\n" 

        # Get wsums for the selected images with commas
        fitsheader *I-model.fits | grep -i wsum | sed s"/WSCVWSUM=\s*//g" > wsums.txt
    fi

    types=("residual" "model")
    mapfile -t sel < selected-channels.txt; export sel;

    for im in ${types[@]}
        do
        echo -e "Copying selected channels' residuls and models for stacking"
        echo -e "---------------------------------------------------------------\n"
        for chan in ${sel[@]}
            do  
                cp  ../*-$chan-I-$im.fits $IMGS
            done

        images=$(ls -v $IMGS/*-I-$im.fits)
        fitstool.py --stack=$SEL_CUBES/i-$im-cube.fits:FREQ $(echo $images)
        done


    # echo "Rename the convolved images"
    # # Adding || tru so that the error here does not fail the entire program
    # # see: https://stackoverflow.com/questions/11231937/bash-ignoring-error-for-a-particular-command
    # rename.ul -- ".convolved.fits" ".fits" $CONV_CUBES/* || true


    echo -e "\n############################################################"
    echo "Do the SPI fitting"
    echo -e "############################################################\n" 

    echo "Normalize the wsums by the largest values"
    echo -e "---------------------------------------------------------------\n"

    # # Doing this with a quick python script because,
    # # well I can :) and store in this variable
    wsums=$(python -c "import numpy as np; wsums = np.loadtxt('wsums.txt'); wsums = np.round(wsums/wsums.max(), 4); print(*wsums)")


    mapfile -t freqs < frequencies.txt

    # cw - channel weights, th-rms threshold factor, 
    # acr - add conv residuals, bm -  beam model
    spimple-spifit -model $SEL_CUBES/i-model-cube.fits \
        -residual $SEL_CUBES/i-residual-cube.fits -o $SPIS/spi-map -th 10 \
        -nthreads 32 -pb-min 0.15 -cw $wsums -acr -bm JimBeam -band l \
        --products aeikb -cf ${freqs[@]}


    return 0
}




generateRmMap(){
    # arg1: str
    #   name: prefix
    #   Prefix prepended to the outputs from this script. 
    #   May include the dump dir. Default is 'initial'.
    
    # local prefix=${1:-"$PRODS/initial"}
    # local mask=${2:-"$mask_dir/true_mask.fits"}
    local args=${1:-"-h"}

    # print function information
    echo "Running function: $FUNCNAME($@)"

    echo -e "\n############################################################"
    echo "Do some RM maps, fpol maps and other maps for each pixel"
    echo -e """Using my mask here, Don't know where yours is but if this step
        \r fails, check on that"""
    echo "Default maximum depth and number of iterations same as that of previous"
    echo -e "############################################################\n"
    sc-rmmap $(echo $args)

    return 0
}


makeMasks(){
    # Automatically create masks
    # arg1: str
    #     Arguments to be passed to houdini

    local args=${1:="-h"}


    echo -e "\n############################################################"
    echo "Running function: $FUNCNAME($args)"
    echo -e "\n############################################################"

    # Make mask for my source
    sc-houdini $(echo $args)
    
    return 0
}



runScrappy(){
    # Run the scrappy package
    # -----------------------
    # arg1: int
    #     name: thresh
    #     Threshold to be used. Default 800
    # arg2:
    #     name: regsize
    #     Region size to be used. Default 3
    # arg3: str
    #     name: scout
    #     Name of output directory for use. Default scrapy-out
    # arg5: bool
    #     name: bokeh
    #     Whether or not to Generate bokeh plots


    # for scrappy
    # I change region size from 5pix to 3 pixels.
    # Set minimum noise * threshold above which to generate regions

    local thresh=${1:-10}
    local regsize=${2:-3}
    local scout=${3:-"$PRODS/scrappy-out"}
    local bokeh=${4:-false}
    local others=${5:-""}

    # print function information
    echo "Running function: $FUNCNAME($@)";
    
    sc-los -rs $regsize --threshold $thresh -odir $scout \
        -ref-image i-mfs.fits -nri i-mfs.fits \
        $(echo $others)

    echo -e "\n############################################################"
    echo "Edit the reg file in a way that it can be loaded into CARTA"
    echo -e "############################################################\n"

    cp $scout/regions/*valid.reg $scout/regions/beacons.reg
    sed -i "s/text=.*//g" $scout/regions/beacons.reg


    echo -e "\n############################################################"
    echo -e """Perfrom RM synthesis for various lines of sight generated from \
        \r previous step and plot the output"""
    echo -e """For pictor I set the maximum depth to 400, depth step 1 looks \
        \r smoothest, niter from 500-1000 looks similar"""
    echo -e "############################################################\n"
    sc-losrm -id $scout/los-data \
        -od $scout/los-rm-data -md 400 --depth-step 1

    if $bokeh
    then
        echo -e "\n############################################################"
        echo "Generate interactive plots for the various LoS"
        echo -e "############################################################\n"
        sc-bokehplot -id $scout/los-data \
          --yaml bk_plots.yml -od $scout/bokeh-plots
    fi

    return 0
}


imageXMask(){

    # Multiply  some image by some mask
    # ima: str
    #     Name of the image
    # mask: str
    #     Mask name
    # output: str
    #     Name of the resulting output image


    local ima=$1
    local mask=$2
    local output=$3

    echo -e "\n############################################################"
    echo -e "Running $FUNCNAME"
    echo -e "############################################################\n"

    #  we premultiply the mask; cubes become problematic otherwise
    fitstool.py --prod $mask $ima -o $output

    return 0;

}

requiredSetup(){
    installRequiredSoftware
    initialiseEnvVarsForThisScript

    return 0
}


# compositeMaskTest(){
#     # ---------------------------------------------------------------------
#     # creating composite mask
#     makeMasks "
#         $PRODS/initial-clean-fdf-SNR.fits
#         -o $mask_dir/fp_snr_mask.fits
#         -above 10
#         "
    
#     # imageXMask(inimage, mask, outimage)
#     imageXMask $PRODS/initial-FPOL-at-max-lpol.fits \
#         $mask_dir/fp_snr_mask.fits \
#         $PRODS/masked-initial-FPOL-at-max-lpol.fits

#     makeMasks "
#         $PRODS/masked-initial-FPOL-at-max-lpol.fits
#         -o $mask_dir/fp_mask.fits
#         -above 0 -below 0.85
#         "
#     # ---------------------------------------------------------------------

#     return 0
# }



main(){
    requiredSetup

    makeDirs
    selectGoodChannels
    stackAllImages

    selectedChannels_Stack
    # selectedChannels_CubeConvolve
    selectedChannels_SinglesConvolve

    copyMfsImages
    
    makeMasks "
        i-mfs.fits
        -o $mask_dir/true_mask.fits
        -above 4e-3
        -rb pica_box_region.reg
        "


    # ----------------------------------------------------------------------
    # select whether to use the convolved images/cubes or not
    CONV=false

    if $CONV
    then
        CUBES=$CONV_CUBES
        IMAGES=$CONVIM
    else
        CUBES=$SEL_CUBES
        IMAGES=$IMGS
    fi
    # ----------------------------------------------------------------------


    # pass all the arguments between the quotes to this function
    generateRmMap "
        -q $CUBES/q-image-cube.fits
        -u $CUBES/u-image-cube.fits
        -i $CUBES/i-image-cube.fits -ncore 45 -md 200
        -o $PRODS/initial -mask $mask_dir/true_mask.fits -f frequencies.txt "

    
    # create a mask for only where fpol SNR > 10
    makeMasks "
        -o $mask_dir/fpol-snrmask.fits 
        -above 10 $PRODS/initial-clean-fdf-SNR.fits"



    # signature
    # runScrappy(thresh, region_size, output_dir, boke_plots?, other_args)
    runScrappy 50 3 $PRODS/scrap-outputs false \
        "-m $mask_dir/fpol-snrmask.fits 
        -nrf pica_noise_region.reg 
        -idir $IMAGES"


    sc-depol $CUBES/i-image-cube.fits \
        $CUBES/q-image-cube.fits \
        $CUBES/u-image-cube.fits \
        $PRODS/scrap-outputs/los-data/reg_1.npz \
        $mask_dir/true_mask.fits

    # generateSpiMap

    return 0
}


LOGFILE="showrunner.log"

if [[ $1 = '-h' ]]
then
    welcome | tee -a $LOGFILE
else
    # Run this main function
    main | tee -a $LOGFILE
fi
