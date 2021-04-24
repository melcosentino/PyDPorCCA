"""
Generates settings file for DPorCCA in txt format.
Parameters:
    Detector:
    - Hydrophone
    - Serial number
    - 
    Classifier:
    - LQ: float
    - HQ: float

    Click Trains:
    - Latitude: in decimal degrees
    - Longitude: in decimal degrees

"""


def generate_set_file(function):
    if function == "detector_porcc":
        f = open("SettingsDPorCCA.txt", "w")
        f.write("Filter: ")
        f.close()
    elif function == "detector_only":
        f = open("SettingsDPorCCA.txt", "w")
        f.write("Woops! I have deleted the content!")
        f.close()
    elif function == "click_trains":
        f = open("SettingsDPorCCA.txt", "w")
        f.write("Latitude: , Longitude")
        f.close()
