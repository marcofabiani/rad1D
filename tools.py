import numpy as np
import re
import pandas as pd

mw_h2o = 18.015 # molecular weight of H2O
mw_co2 = 44.01 # molecular weight of CO2
sigma = 5.67e-8 # Stefan-Boltzmann constant

# computes weighted-sum-of-gray-gases emissivity 
# T = temperature in K
# p = pressure in bar
# L = characteristic length in m
# xH2O = mole fraction of H2O
# xCO2 = mole fraction of CO2
def wsgg(T,p,L,xH2O,xCO2):
  
  Tref = 2300.0 
  if (T<1000.): 
    T = 1000.
  
  MR = xH2O/(xCO2+1e-5)
  MRvec = np.zeros(11)
  MRvec[0:11] = [0.125,0.25,0.5,0.75,1.,2.,2.5,3,4.,6.,8.]

  idx = (np.abs(MRvec - MR)).argmin()

  k = np.zeros(5)
  cij = np.zeros((4,4))
  
  if (idx==0) :
    k = [  0.201355,  83.681211,   2.177752,   0.018921]
    cij[0,:] = [  0.338670,   0.126308,  -0.373209,   0.123367 ] 
    cij[1,:] = [  0.198925,  -0.301228,   0.165138,  -0.031874] 
    cij[2,:] = [  0.582898,  -0.935571,   0.565397,  -0.120346] 
    cij[3,:] = [ -0.532591,   2.300377,  -1.761992,   0.402114]
  elif (idx==1):
    k = [   0.176943,  68.107015,   1.927536,   0.019695] 
    cij[0,:] = [  0.013091,   0.948018,  -0.939157,   0.244841] 
    cij[1,:] = [  0.252767,  -0.430343,   0.267167,  -0.058113] 
    cij[2,:] = [  0.728348,  -1.080919,   0.581898,  -0.108590] 
    cij[3,:] = [ -0.396877,   1.790107,  -1.318230,   0.291051] 
  elif (idx==2):
    k = [   0.167276,   1.952047,  52.441380,   0.018473] 
    cij[0,:] = [ -0.206884,   1.486154,  -1.253509,   0.298121] 
    cij[1,:] = [  0.725780,  -0.884108,   0.345353,  -0.036283] 
    cij[2,:] = [  0.335381,  -0.627902,   0.419616,  -0.096346] 
    cij[3,:] = [ -0.236388,   1.262742,  -0.904787,   0.195389] 
  elif (idx==3):
    k = [  0.169308,  42.244993,   1.969425,   0.017736] 
    cij[0,:] = [ -0.277568,   1.632683,  -1.308007,   0.298863] 
    cij[1,:] = [  0.401490,  -0.778331,   0.531477,  -0.123635] 
    cij[2,:] = [  0.670017,  -0.660894,   0.137143,   0.020843] 
    cij[3,:] = [ -0.159840,   1.029515,  -0.728327,   0.155358] 
  elif (idx==4):
    k = [   0.174844,  36.680098,   1.995313,   0.017628] 
    cij[0,:] = [ -0.302002,   1.667088,  -1.300099,   0.289894] 
    cij[1,:] = [  0.450271,  -0.887988,   0.611855,  -0.142996] 
    cij[2,:] = [  0.617410,  -0.486243,  -0.015639,   0.061352] 
    cij[3,:] = [ -0.120589,   0.914228,  -0.640886,   0.135256] 
  elif (idx==5):
    k = [   0.186016,  25.467792,   1.950971,   0.017696] 
    cij[0,:] = [ -0.320818,   1.637206,  -1.207716,   0.254019] 
    cij[1,:] = [  0.584577,  -1.174975,   0.814065,  -0.190224] 
    cij[2,:] = [  0.458295,  -0.017773,  -0.401982,   0.160087] 
    cij[3,:] = [ -0.056075,   0.724775,  -0.497380,   0.102112] 
  elif (idx==6):
    k = [  0.189185,  23.112324,   1.928168,   0.017814] 
    cij[0,:] = [ -0.318111,   1.608209,  -1.169021,   0.241602] 
    cij[1,:] = [  0.625255,  -1.259805,   0.872198,  -0.203453] 
    cij[2,:] = [  0.406785,   0.124317,  -0.514596,   0.188062] 
    cij[3,:] = [ -0.043536,   0.687196,  -0.468620,   0.095378] 
  elif (idx==7):
    k = [   0.191435,  21.511787,   1.908193,   0.017914] 
    cij[0,:] = [ -0.314168,   1.581880,  -1.137337,   0.231845] 
    cij[1,:] = [  0.656536,  -1.324673,   0.916222,  -0.213366] 
    cij[2,:] = [  0.366303,   0.234028,  -0.600378,   0.209141] 
    cij[3,:] = [ -0.035419,   0.662399,  -0.449544,   0.090890] 
  elif (idx==8):
    k = [   0.194536,  19.480684,   1.878634,   0.018078] 
    cij[0,:] = [ -0.305586,   1.538113,  -1.089016,   0.217523] 
    cij[1,:] = [  0.701315,  -1.417862,   0.979208,  -0.227454] 
    cij[2,:] = [  0.306826,   0.392967,  -0.723220,   0.239031] 
    cij[3,:] = [ -0.026015,   0.632872,  -0.426646,   0.085466] 
  elif (idx==9):
    k = [  0.198014,  17.336022,   1.842568,   0.018294] 
    cij[0,:] = [ -0.290310,   1.475621,  -1.025942,   0.199634] 
    cij[1,:] = [  0.755432,  -1.533041,   1.057727,  -0.245063] 
    cij[2,:] = [  0.232303,   0.590153,  -0.874267,   0.275481] 
    cij[3,:] = [ -0.018249,   0.606954,  -0.406312,   0.080618] 
  elif (idx==10):
    k = [   0.199856,  16.157754,   1.820514,   0.018416] 
    cij[0,:] = [ -0.278339,   1.432731,  -0.985239,   0.188458] 
    cij[1,:] = [  0.788192,  -1.605694,   1.108504,  -0.256651] 
    cij[2,:] = [  0.185639,   0.713984,  -0.969253,   0.298410] 
    cij[3,:] = [ -0.015392,   0.596299,  -0.397852,   0.078594] 
  
  epsg = 0
  ig = 0
  a = np.zeros(len(k))

  for kg in k:
    j=0
    for cc in cij[ig,:]:
      a[ig]=a[ig]+cij[ig,j]*(T/Tref)**(j)
    
      j = j+1  
      

    kk = k[ig]*(xH2O+xCO2)*p
   
    epsg=epsg+a[ig]*(1-np.exp(-kk*L))

    ig = ig+1 
    


  return epsg

# routines to extract data from CEA output
_num = re.compile(r'^[\+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][\+\-]?\d+)?$')

def _is_number(tok: str) -> bool:
    return bool(_num.match(tok))

def _merge_broken_sci(tokens):
    merged = []
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if (i + 1 < len(tokens)
            and re.match(r'^[\+\-]?(?:\d+(?:\.\d*)?|\.\d+)$', t)
            and tokens[i+1] in ("0", "+0", "-0")):
            merged.append(t + "e" + tokens[i+1])
            i += 2
        else:
            merged.append(t)
            i += 1
    return merged

def _parse_three_col_block(lines, header_regex, stop_regexes):
    start = None
    for i, line in enumerate(lines):
        if re.search(header_regex, line):
            start = i
            break
    if start is None:
        return {}

    records = {}
    i = start + 1
    while i < len(lines) and not any(re.search(p, lines[i]) for p in stop_regexes):
        line = lines[i]
        if not line.strip():  # <-- previously skipped blank lines
            i += 1
            continue
        parts = line.split()

        # Build label until numbers begin
        j = 0
        label_tokens = []
        while j < len(parts):
            t = parts[j]
            broken_pair = (j + 1 < len(parts)
                           and re.match(r'^[\+\-]?(?:\d+(?:\.\d*)?|\.\d+)$', t)
                           and parts[j+1] in ("0", "+0", "-0"))
            if _is_number(t) or broken_pair:
                break
            label_tokens.append(t)
            j += 1
        label = " ".join(label_tokens) if label_tokens else parts[0]
        label = label.rstrip(",")

        # parse numbers
        num_tokens = _merge_broken_sci(parts[j:])
        nums = []
        for t in num_tokens:
            try:
                nums.append(float(t))
            except ValueError:
                pass
        if len(nums) >= 3:
            records[label] = nums[:3]

        i += 1

    return records

def extract_chamber_throat_exit_and_molefractions(output_str, normalize_keys=True):
    lines = output_str.splitlines()

    perf = _parse_three_col_block(
        lines,
        header_regex=r'CHAMBER\s+THROAT\s+EXIT',
        stop_regexes=[
            r'^\s*TRANSPORT PROPERTIES',
            r'^\s*WITH EQUILIBRIUM REACTIONS',
            r'^\s*WITH FROZEN REACTIONS',
            r'^\s*PERFORMANCE PARAMETERS',
            r'^\s*MOLE FRACTIONS',
        ],
    )

    mole = _parse_three_col_block(
        lines,
        header_regex=r'^\s*MOLE FRACTIONS\s*$',
        stop_regexes=[
            r'^\s*\* THERMODYNAMIC PROPERTIES',
            r'^\s*PRODUCTS WHICH WERE CONSIDERED',
            r'^\s*NOTE\.',
        ],
    )

    all_records = {**perf, **mole}

    if normalize_keys:
        cleaned = {}
        for k, v in all_records.items():
            kk = k.lstrip('*')
            kk = re.sub(r'\s+', ' ', kk)
            cleaned[kk] = v
        all_records = cleaned

    df = pd.DataFrame(all_records, index=["CHAMBER", "THROAT", "EXIT"]).T
    return df