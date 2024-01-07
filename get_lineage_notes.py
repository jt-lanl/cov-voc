import urllib.request
LINK="https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt"
RAWLINK="https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
try:
   with urllib.request.urlopen(RAWLINK) as f:
      print(f.read().decode('utf-8'))
except urllib.error.URLError as e:
   print(e.reason)
