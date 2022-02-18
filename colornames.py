## Here are the so-called "x11" color names
## You can see the colors at https://en.wikipedia.org/wiki/X11_color_names

import re

COLORNAMES='''#F0F8FF AliceBlue
#FAEBD7 AntiqueWhite
#00FFFF Aqua
#7FFFD4 Aquamarine
#F0FFFF Azure
#F5F5DC Beige
#FFE4C4 Bisque
#000000 Black
#FFEBCD BlanchedAlmond
#0000FF Blue
#8A2BE2 BlueViolet
#A52A2A Brown
#DEB887 BurlyWood
#5F9EA0 CadetBlue
#7FFF00 Chartreuse
#D2691E Chocolate
#FF7F50 Coral
#6495ED CornflowerBlue
#FFF8DC Cornsilk
#DC143C Crimson
#00FFFF Cyan
#00008B DarkBlue
#008B8B DarkCyan
#B8860B DarkGoldenrod
#A9A9A9 DarkGray
#006400 DarkGreen
#BDB76B DarkKhaki
#8B008B DarkMagenta
#556B2F DarkOliveGreen
#FF8C00 DarkOrange
#9932CC DarkOrchid
#8B0000 DarkRed
#E9967A DarkSalmon
#8FBC8F DarkSeaGreen
#483D8B DarkSlateBlue
#2F4F4F DarkSlateGray
#00CED1 DarkTurquoise
#9400D3 DarkViolet
#FF1493 DeepPink
#00BFFF DeepSkyBlue
#696969 DimGray
#1E90FF DodgerBlue
#B22222 FireBrick
#FFFAF0 FloralWhite
#228B22 ForestGreen
#FF00FF Fuchsia
#DCDCDC Gainsboro
#F8F8FF GhostWhite
#FFD700 Gold
#DAA520 Goldenrod
#808080 Gray
#008000 Green
#ADFF2F GreenYellow
#F0FFF0 Honeydew
#FF69B4 HotPink
#CD5C5C IndianRed
#4B0082 Indigo
#FFFFF0 Ivory
#F0E68C Khaki
#E6E6FA Lavender
#FFF0F5 LavenderBlush
#7CFC00 LawnGreen
#FFFACD LemonChiffon
#ADD8E6 LightBlue
#F08080 LightCoral
#E0FFFF LightCyan
#FAFAD2 LightGoldenrodYellow
#90EE90 LightGreen
#D3D3D3 LightGrey
#FFB6C1 LightPink
#FFA07A LightSalmon
#20B2AA LightSeaGreen
#87CEFA LightSkyBlue
#778899 LightSlateGray
#B0C4DE LightSteelBlue
#FFFFE0 LightYellow
#00FF00 Lime
#32CD32 LimeGreen
#FAF0E6 Linen
#FF00FF Magenta
#800000 Maroon
#66CDAA MediumAquamarine
#0000CD MediumBlue
#BA55D3 MediumOrchid
#9370DB MediumPurple
#3CB371 MediumSeaGreen
#7B68EE MediumSlateBlue
#00FA9A MediumSpringGreen
#48D1CC MediumTurquoise
#C71585 MediumVioletRed
#191970 MidnightBlue
#F5FFFA MintCream
#FFE4E1 MistyRose
#FFE4B5 Moccasin
#FFDEAD NavajoWhite
#000080 Navy
#FDF5E6 OldLace
#808000 Olive
#6B8E23 OliveDrab
#FFA500 Orange
#FF4500 OrangeRed
#DA70D6 Orchid
#EEE8AA PaleGoldenrod
#98FB98 PaleGreen
#AFEEEE PaleTurquoise
#DB7093 PaleVioletRed
#FFEFD5 PapayaWhip
#FFDAB9 PeachPuff
#CD853F Peru
#FFC0CB Pink
#DDA0DD Plum
#B0E0E6 PowderBlue
#800080 Purple
#FF0000 Red
#BC8F8F RosyBrown
#4169E1 RoyalBlue
#8B4513 SaddleBrown
#FA8072 Salmon
#F4A460 SandyBrown
#2E8B57 SeaGreen
#FFF5EE Seashell
#A0522D Sienna
#C0C0C0 Silver
#87CEEB SkyBlue
#6A5ACD SlateBlue
#708090 SlateGray
#FFFAFA Snow
#00FF7F SpringGreen
#4682B4 SteelBlue
#D2B48C Tan
#008080 Teal
#D8BFD8 Thistle
#FF6347 Tomato
#40E0D0 Turquoise
#EE82EE Violet
#F5DEB3 Wheat
#FFFFFF White
#F5F5F5 WhiteSmoke
#FFFF00 Yellow
#9ACD32 YellowGreen'''

import random

## upon import, create ColorTable dict
ColorTable = dict()
for hn in COLORNAMES.splitlines():
    hexname,commonname = hn.strip().split()
    ColorTable[commonname.lower()] = hexname

def tohex(commonname):
    '''convert commonname (eg, YellowGreen) to hex string;
    if common name is given as a hexstring already, then return the hexstring.
    valid input hexstrings of form #9ACD32 or 9ACD32
    output hexstring will always be preceeded by #
    '''
    name = commonname.lower()
    try:
        return ColorTable[name]
    except KeyError:
        if re.match(r'\#[0-9a-f]{6}$',name):
            return name
        elif re.match(r'[0-9a-f]{6}$',name):
            return "#"+name
        else:
            raise KeyError('Invalid color: {name}')

def torgb(commonname):
    hex = tohex(commonname)
    rgb = [int("0x" + c, 16) for c in (hex[1:3],hex[3:5],hex[5:7])]
    return rgb

def lighter(commonname):
    '''return a lighter color than the input color'''
    r,g,b = torgb(commonname)
    r = (r + 255) // 2
    g = (g + 255) // 2
    b = (b + 255) // 2
    return '#' + "".join([format(v,"02x") for v in (r,g,b)])

def darker(commonname):
    '''return a darker color than the input color'''
    r,g,b = torgb(commonname)
    r,g,b = [v//2 for v in (r,g,b)] 
    return '#' + "".join([format(v,"02x") for v in (r,g,b)])
    
def random_hex():
    return random.choice(list(ColorTable.values()))

def random_color():
    return random.choice(list(ColorTable.keys()))

if __name__ == "__main__":

    for t in ColorTable:
        print(f"{t:20s} {ColorTable[t]}")

    for t in ['Yellow', '06A6f6', '#06a6F6']:
        print(f"{t:20s} {tohex(t)}")

    name = 'YellowOrange'
    name = '#06a6F77'
    try:
        hex=tohex(name)
        print(f"{name:20s} {hex}")
    except Exception as e:
        print(f"Failed, as it should have with name={name}:",e)
        
        
        
