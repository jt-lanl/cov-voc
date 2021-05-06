lightgreen  = '#77b300'
mediumgreen = '#558000'
darkgreen   = '#334d00'
lightred    = '#bf00ff'
darkred     = '#8600b3'
black       = '#000000'

_colordict = {
    14: lightgreen,
    15: lightgreen,
    16: lightgreen,
    17: lightgreen,
    18: lightgreen,
    19: lightgreen,
    20: lightgreen,
    140: mediumgreen,
    141: mediumgreen,
    142: mediumgreen,
    143: mediumgreen,
    144: mediumgreen,
    145: mediumgreen,
    146: mediumgreen,
    147: mediumgreen,
    148: mediumgreen,
    149: mediumgreen,
    150: mediumgreen,
    151: mediumgreen,
    152: mediumgreen,
    153: mediumgreen,
    154: mediumgreen,
    155: mediumgreen,
    156: mediumgreen,
    157: mediumgreen,
    158: mediumgreen,
    242: darkgreen,
    243: darkgreen,
    244: darkgreen,
    245: darkgreen,
    246: darkgreen,
    247: darkgreen,
    248: darkgreen,
    249: darkgreen,
    250: darkgreen,
    251: darkgreen,
    252: darkgreen,
    253: darkgreen,
    254: darkgreen,
    255: darkgreen,
    256: darkgreen,
    257: darkgreen,
    258: darkgreen,
    259: darkgreen,
    260: darkgreen,
    261: darkgreen,
    262: darkgreen,
    263: darkgreen,
    264: darkgreen,
    330: lightred,
    331: lightred,
    332: lightred,
    333: lightred,
    334: lightred,
    335: lightred,
    336: lightred,
    337: lightred,
    338: lightred,
    339: lightred,
    340: lightred,
    341: lightred,
    342: lightred,
    343: lightred,
    344: lightred,
    345: lightred,
    346: lightred,
    347: lightred,
    348: lightred,
    349: lightred,
    350: lightred,
    351: lightred,
    352: lightred,
    353: lightred,
    354: lightred,
    355: lightred,
    356: lightred,
    357: lightred,
    358: lightred,
    359: lightred,
    360: lightred,
    361: lightred,
    362: lightred,
    363: lightred,
    364: lightred,
    365: lightred,
    366: lightred,
    367: lightred,
    368: lightred,
    369: lightred,
    370: lightred,
    371: lightred,
    372: lightred,
    373: lightred,
    374: lightred,
    375: lightred,
    376: lightred,
    377: lightred,
    378: lightred,
    379: lightred,
    380: lightred,
    381: lightred,
    382: lightred,
    383: lightred,
    384: lightred,
    385: lightred,
    386: lightred,
    387: lightred,
    388: lightred,
    389: lightred,
    390: lightred,
    391: lightred,
    392: lightred,
    393: lightred,
    394: lightred,
    395: lightred,
    396: lightred,
    397: lightred,
    398: lightred,
    399: lightred,
    400: lightred,
    401: lightred,
    402: lightred,
    403: lightred,
    404: lightred,
    405: lightred,
    406: lightred,
    407: lightred,
    408: lightred,
    409: lightred,
    410: lightred,
    411: lightred,
    412: lightred,
    413: lightred,
    414: lightred,
    415: lightred,
    416: lightred,
    417: lightred,
    418: lightred,
    419: lightred,
    420: lightred,
    421: lightred,
    422: lightred,
    423: lightred,
    424: lightred,
    425: lightred,
    426: lightred,
    427: lightred,
    428: lightred,
    429: lightred,
    430: lightred,
    431: lightred,
    432: lightred,
    433: lightred,
    434: lightred,
    435: lightred,
    436: lightred,
    437: lightred,
    438: darkred,
    439: darkred,
    440: darkred,
    441: darkred,
    442: darkred,
    443: darkred,
    444: darkred,
    445: darkred,
    446: darkred,
    447: darkred,
    448: darkred,
    449: darkred,
    450: darkred,
    451: darkred,
    452: darkred,
    453: darkred,
    454: darkred,
    455: darkred,
    456: darkred,
    457: darkred,
    458: darkred,
    459: darkred,
    460: darkred,
    461: darkred,
    462: darkred,
    463: darkred,
    464: darkred,
    465: darkred,
    466: darkred,
    467: darkred,
    468: darkred,
    469: darkred,
    470: darkred,
    471: darkred,
    472: darkred,
    473: darkred,
    474: darkred,
    475: darkred,
    476: darkred,
    477: darkred,
    478: darkred,
    479: darkred,
    480: darkred,
    481: darkred,
    482: darkred,
    483: darkred,
    484: darkred,
    485: darkred,
    486: darkred,
    487: darkred,
    488: darkred,
    489: darkred,
    490: darkred,
    491: darkred,
    492: darkred,
    493: darkred,
    494: darkred,
    495: darkred,
    496: darkred,
    497: darkred,
    498: darkred,
    499: darkred,
    500: darkred,
    501: darkred,
    502: darkred,
    503: darkred,
    504: darkred,
    505: darkred,
    506: darkred,
    507: lightred,
    508: lightred,
    509: lightred,
    510: lightred,
    511: lightred,
    512: lightred,
    513: lightred,
    514: lightred,
    515: lightred,
    516: lightred,
    517: lightred,
    518: lightred,
    519: lightred,
    520: lightred,
    521: lightred,
}

_boldlist = [
    18,
    144,
    149,
    242,
    243,
    244,
    246,
    253,
    417,
    441,
    444,
    446,
    450,
    452,
    484,
]

def getcolor(site):
    return _colordict[site] if site in _colordict else black

def getweight(site):
    return "bold" if site in _boldlist else "normal"

if __name__ == "__main__":
    ## tests
    for s in [18, 144, 149, 14, 140, 441, 437, 27]:
        print(s,getcolor(s),getweight(s))
        