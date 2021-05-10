from madrigalWeb import madrigalWeb
import re
import os
import pathlib

if __name__ == '__main__':
    """
    List instruments and experiments from Madrigal
    Filter for World Day Experiments
    Optionally Download Experiment Data
    """

    # Choose which radar to use
    RISRN = False
    RISRC = True

    DOWNLOAD_HDF5 = True  # To prevent downloading when you just want to browse instruments and experiments

    Arr_offset_to_download = 19  # Browse first to find this

    if RISRN and RISRC:
        raise Exception("You cannot choose both RISR-N and C because their data lives in different places")
    if not RISRN and not RISRC:
        raise Exception("You have to to pick one of RISR-N or C")

    # NOTE: Different instruments are at different URLs
    if RISRN:
        madrigalUrl = "http://isr.sri.com/madrigal/"
    else:
        madrigalUrl = "https://madrigal.phys.ucalgary.ca/"
    testData = madrigalWeb.MadrigalData(madrigalUrl)

    # List all of the instruments
    instList = testData.getAllInstruments()
    inst_code = 0
    radar_nmonic = ""
    for inst in instList:
        #  91: Resolute Bay North IS Radar (http://isr.sri.com/madrigal/)
        if inst.code == 91 and RISRN:
            inst_code = 91
            print(inst)
            radar_nmonic = inst.mnemonic
        #  92: Resolute Bay Canada IS Radar (https://madrigal.phys.ucalgary.ca/)
        if inst.code == 92 and RISRC:
            inst_code = 92
            print(inst)
            radar_nmonic = inst.mnemonic
    if inst_code == 0 or radar_nmonic == "":
        raise Exception("Something went wrong, inst_code and radar_nmonic were not set.")

    # Get a list of world day experiments
    worldDayRegex = re.compile("WorldDay*")
    expList = testData.getExperiments(inst_code,       # instrument code
                                      2016, 1, 1,  # start year, month, and day
                                      0, 0, 0,   # start hour, minute and second
                                      2017, 1, 1,  # end year, month, and day
                                      0, 0, 0)   # end hour, minute and second
    worldDayExpList = [exp for exp in expList if re.match(worldDayRegex, exp.name)]
    print("--\nNumber of World Day Experiments here: " + str(len(worldDayExpList)) + "\n")
    for i, exp in enumerate(worldDayExpList):
        print(exp.name)
        print("Experiment id: " + str(exp.id))
        print("Experiment date: " + str(exp.startyear) + "." + str(exp.startmonth) + "." + str(exp.startday) + " - "
              + str(exp.endyear) + "." + str(exp.endmonth) + "." + str(exp.endday))
        print("Array offset: " + str(i) + "\n")

    # Download the data from one of the experiments
    if DOWNLOAD_HDF5:
        expToDownload = worldDayExpList[Arr_offset_to_download]
        print("--\nWe are going to download the following experiments:")
        print(expToDownload.name)
        print("Experiment id: " + str(expToDownload.id))
        print("Experiment date: " + str(expToDownload.startyear) + "." + str(expToDownload.startmonth) + "." +
              str(expToDownload.startday) + " - " + str(expToDownload.endyear) + "." + str(expToDownload.endmonth) + "."
              + str(expToDownload.endday))
        destination = "data/" + radar_nmonic + "/" + radar_nmonic \
                      + str(expToDownload.startyear) + str(expToDownload.startmonth) + str(expToDownload.startday) \
                      + "/" + radar_nmonic + str(expToDownload.startyear) + str(expToDownload.startmonth) + str(expToDownload.startday)
        fileList = testData.getExperimentFiles(expToDownload.id)
        i = 1

        # Ensure out directory
        loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
        out_dir = loc_root + "/DataReading/RISR/data/" + radar_nmonic + "/" + radar_nmonic \
                  + str(expToDownload.startyear) + str(expToDownload.startmonth) + str(expToDownload.startday)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Actually download
        for file in fileList:
            print("Downloading: " + file.name + " as " + destination + "." + str(i) + ".h5")
            testData.downloadFile(file.name, destination + "." + str(i) + ".h5", "Michael Luciuk", "mrl280@usask.ca",
                                  "University of Saskatchewan", format="hdf5")
            i = i + 1
