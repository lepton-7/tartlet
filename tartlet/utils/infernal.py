import requests

from tart.utils.utils import print

current_ver = 14.9

# As of rfam 14.9
rfam_riboswitch_accessions = [
    "RF00050",
    "RF00059",
    "RF00080",
    "RF00162",
    "RF00167",
    "RF00168",
    "RF00174",
    "RF00234",
    "RF00379",
    "RF00380",
    "RF00442",
    "RF00504",
    "RF00521",
    "RF00522",
    "RF00634",
    "RF01054",
    "RF01055",
    "RF01056",
    "RF01057",
    "RF01482",
    "RF01510",
    "RF01689",
    "RF01704",
    "RF01725",
    "RF01727",
    "RF01734",
    "RF01739",
    "RF01750",
    "RF01767",
    "RF01786",
    "RF01826",
    "RF01831",
    "RF02680",
    "RF02683",
    "RF02885",
    "RF02912",
    "RF03057",
    "RF03058",
    "RF03071",
    "RF03072",
]


def make_riboswitch_cm():
    with open("data/riboswitches.cm", "w") as f:
        for acc in rfam_riboswitch_accessions:
            address = f"https://rfam.org/family/{acc}/cm"

            r = requests.get(address)

            f.write(r.text)
            print(f"Ingested {acc}")


def get_clanin(ver: str = current_ver):
    url = f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{ver}/Rfam.clanin"

    with open(f"data/rfam_{ver}.clanin", "w") as f:
        r = requests.get(url)
        f.write(r.text)
