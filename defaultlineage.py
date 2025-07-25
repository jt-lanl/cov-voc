DefaultLineageTable = [

    (       'Black',         'Unassigned', r'(None|Unassigned|)'),
    (    'DarkGray',          'Ancestral', r'(A)|(A\.2)|(A\.2\.[^5](\..*)?)|(A.[^2](\..*)?)|(A\.1[01234567](\..*)?)|(A\.[^1][0-9](\..*)?)|(B)|(B\.1\.1\.161)|(B\.1\.14)|(B\.1\.260)|(B\.1[0123589])|(B\.2[036789])|(B\.3(\.1)?)|(B.3[0-9])|(B\.4(\.[124567])?)|(B\.4[012345679])|(B\.5([01235678]?))|(B.6(\.[1234568])?)|(B\.6[01])'),
    ('AntiqueWhite', 'Other_Recombinants', r'X.*'),
    (      'Yellow',        'D614G/Other', r'$OTHER'),
    (      'Orange',              'Alpha', r'((B\.1\.1\.7)|(Q))(\..*)?'),
    (    'DarkCyan',               'Beta', r'B\.1\.351(\..*)?'),
    (   'FireBrick',              'Gamma', r'((B\.1\.1\.28)|(P))(\..*)?'),
    (       'Green',            'Epsilon', r'B\.1\.42[97](\..*)?'),
    (   'LightBlue',               'Iota', r'B\.1\.526(\..*)?'),
    (  'BlueViolet',              'Delta', r'((B\.1\.617\.2)|(AY))(\..*)?'),
    (     'HotPink',               'BA.1', r'((BA\.1)|(B[CD]))(\..*)?'),
    ('LightSkyBlue',               'BA.2', r'((BA\.2)|(B[GHJLMNPRSY])|(C[ABHJMV])|(D[DSV])|(E[JP])|(F[JKRSV])|(G[PQ])|(J[LNPQV])|(K[GPQRSUVWYZ])|(L[ABCDEFGHJKLMNPQRSTUVWYZ])|(M[ABCDEFGHJKLMNPQRSTUVWYZ])|(N[ACDEFGHJKLMNPQRSTUVWYZ])|(P[ABCDEFGHJKLMPRSVWYZ])|(Q[ABC]))(\..*)?'),
    (         'Red',          'BA.2.12.1', r'((BA\.2\.12\.1)|(BG))(\..*)?'),
    (   'Chocolate',               'BA.4', r'((BA\.4)|(CS)|(DC))(\..*)?'),
    ('MidnightBlue',               'BA.5', r'((BA\.5)|(B[EFKQTUVWZ])|(C[CDEFGKLNPQRTUWYZ])|(D[ABEFGHJKLMNPQRTUWYZ])|(E[ABCDEFHNQRSTVWYZ])|(F[ABCFMNQ])|(JH))(\..*)?'),
    (       'Khaki',           'BQ.1.1.1', r'((BQ\.1\.1\.1)|(CZ))(\..*)?'),
    ( 'ForestGreen',            'BA.2.75', r'((BA\.2\.75)|(B[LMNRY])|(C[ABHJV])|(D[SV])|(E[JP])|(F[JKRS])|(G[PQ])|(J[LPV])|(KG))(\..*)?'),
    (     'DarkRed',                'XBB', r'((XBB)|(E[GKLMU])|(F[DEGHLPTUWYZ])|(G[ABCDEFGHJKMNRSUVWYZ])|(H[ABCDEFGHJKLMNPQRSTUVYZ])|(J[ABCDEFGJKMRSUWYZ])|(K[ABCEFHJKLMNT]))(\..*)?'),
    (      'Salmon',            'XBB.1.9', r'((XBB\.1\.9)|(EG)|(FL)|(GD)|(H[KNV])|(J[GJR])|(K[BCFL]))(\..*)?'),
    (     'Crimson',            'XBB.1.5', r'((XBB\.1\.5)|(E[KLMU])|(F[DGHTZ])|(G[BCFGKNRUV])|(H[ACDJMPQRSTYZ])|(J[BDKZ])|(K[AKMN]))(\..*)?'),
    (   'LimeGreen',             'CH.1.1', r'((CH\.1\.1)|(DV)|(F[JKS])|(G[PQ])|(J[LPV])|(KG))(\..*)?'),
    (      'Purple',             'EG.5.1', r'((EG\.5\.1)|(H[KV])|(J[GJR])|(K[BL]))(\..*)?'),
    (     'Thistle',               'HK.3', r'HK\.3(\..*)?'),
    (      'Violet',               'HV.1', r'((HV\.1)|(KL))(\..*)?'),
    (  'Aquamarine',            'XBB.2.3', r'((XBB\.2\.3)|(G[EJMSZ])|(H[GH])|(J[AESUY])|(K[HT]))(\..*)?'),
    (  'DodgerBlue',            'BA.2.86', r'((BA\.2\.86)|(J[NQ])|(K[PQRSUVWYZ])|(L[ABCDEFGHJKLMNPQRSTUVWYZ])|(M[ABCDEFGHJKLMNPQRSTUVWYZ])|(N[ACDEFGHJKLMNPQRSTUVWYZ])|(P[ABCDEFGHJKLMPRSVWYZ])|(Q[ABC]))(\..*)?'),
    (        'Cyan',               'JN.1', r'((JN\.1)|(K[PQRSUVWYZ])|(L[ABCDEFGHJKLMNPQRSTUWYZ])|(M[ABCDEFGHJKLMNQRSTUVWYZ])|(N[ACDEFGHJKLMNPQRSTUVWYZ])|(P[ABCDEFGHJKLMPRSVWYZ])|(Q[ABC]))(\..*)?'),
    (       'Brown',           'KP.1.1.3', r'((KP\.1\.1\.3)|(LP)|(N[WY])|(P[DFPR]))(\..*)?'),
    (   'RoyalBlue',               'KP.2', r'((KP\.2)|(MW)|(N[GHKMN]))(\..*)?'),
    (    'DarkBlue',             'KP.2.3', r'((KP\.2\.3)|(MW)|(N[GKMN]))(\..*)?'),
    (        'Plum',               'KP.3', r'((KP\.3)|(LW)|(M[CKLMRY])|(N[PQRV])|(P[ABEGHJM]))(\..*)?'),
    (     'Magenta',           'KP.3.1.1', r'((KP\.3\.1\.1)|(MC)|(P[AEHJ]))(\..*)?'),
    (   'Turquoise',                'XEC', r'((XEC)|(P[TU]))(\..*)?'),
    (   'PaleGreen',             'XEC.25', r'XEC\.25(\..*)?'),
    (       'Khaki',             'LP.8.1', r'((LP\.8\.1)|(N[WY])|(P[DFPR]))(\..*)?'),
    (       'Coral',           'LP.8.1.1', r'((LP\.8\.1\.1)|(NY))(\..*)?'),
    (   'FireBrick',               'LF.7', r'((LF\.7)|(NT)|(P[CLVWYZ])|(Q[ABC]))(\..*)?'),
    ( 'ForestGreen',                'XDV', r'((XDV)|(NB)|(P[NQ]))(\..*)?'),
    (   'LawnGreen',           'NB.1.8.1', r'((NB\.1\.8\.1)|(PQ))(\..*)?'),
    (  'DodgerBlue',                'XFG', r'XFG(\..*)?'),
    (       'Black',                'XFC', r'XFC(\..*)?'),

]
