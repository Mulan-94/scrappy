from colorama import init, Fore, Style

init(autoreset=True)

def avail_color_functions():
    colours = "BLUE CYAN GREEN MAGENTA RED WHITE YELLOW".split()

    def format_to_colour(colour, bold=False):
        if bold:
            return lambda string: f"{Fore.__getattribute__(colour)}{Style.BRIGHT}{string}{Style.RESET_ALL}"
        else:
            return lambda string: f"{Fore.__getattribute__(colour)}{string}{Style.RESET_ALL}"

    colour_funcs = {colour.lower(): format_to_colour(colour) for colour in colours}
    colour_funcs.update({
        f"bold_{colour.lower()}": format_to_colour(colour, bold=True) for colour in colours})
  
    return colour_funcs



