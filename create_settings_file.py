"""
Generates settings file for DPorCCA in txt format.
Parameters:
    Detector:
    - Hydrophone
    - Serial number
    - 
    Classifier:
    - LQ
    - HQ
    Click Trains:
    - Latitude
    - Longitude

"""



def generate_set_file(function):
    if function == "Detector":
        f = open("SettingsDPorCCA.txt", "w")
        f.write("Woops! I have deleted the content!")
        f.close()
