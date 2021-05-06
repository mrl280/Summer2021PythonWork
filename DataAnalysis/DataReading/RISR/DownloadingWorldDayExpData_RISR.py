from madrigalWeb import madrigalWeb
import re

if __name__ == '__main__':
    """
    List instruments and experiments from Madrigal
    Filter for World Day Experiments
    Optionally Download Experiment Data
    """
    # To prevent downloading when you just want to browse instruments and experiments
    DOWNLOAD_HDF5 = False

    madrigalUrl = "http://isr.sri.com/madrigal/"  # NOTE: Different instruments are at different URLs
    testData = madrigalWeb.MadrigalData(madrigalUrl)

    # List all of the instruments
    instList = testData.getAllInstruments()
    for inst in instList:
        #  91: Resolute Bay North IS Radar (http://isr.sri.com/madrigal/)
        if inst.code == 91:
            print(inst)
        #  92: Resolute Bay Canada IS Radar (https://madrigal.phys.ucalgary.ca/)
        # if inst.code == 92:
        #     print(inst)

    # Get a list of world day experiments
    worldDayRegex = re.compile("WorldDay*")
    expList = testData.getExperiments(91,       # instrument code
                                      2016, 1, 1,  # start year, month, and day
                                      0, 0, 0,   # start hour, minute and second
                                      2018, 1, 1,  # end year, month, and day
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
        expToDownload = worldDayExpList[6]
        print("--\nWe are going to download the following experiments:")
        print(expToDownload.name)
        print("Experiment id: " + str(expToDownload.id))
        print("Experiment date: " + str(expToDownload.startyear) + "." + str(expToDownload.startmonth) + "." +
              str(expToDownload.startday) + " - " + str(expToDownload.endyear) + "." + str(expToDownload.endmonth) + "."
              + str(expToDownload.endday))
        destination = "data/ran" + str(expToDownload.startyear) + str(expToDownload.startmonth) + str(expToDownload.startday) \
                      + "/ran" + str(expToDownload.startyear) + str(expToDownload.startmonth) + str(expToDownload.startday)
        fileList = testData.getExperimentFiles(expToDownload.id)
        i = 1
        for file in fileList:
            print("Downloading: " + file.name + " as " + destination + "." + str(i) + ".h5")
            testData.downloadFile(file.name, destination + "." + str(i) + ".h5", "Michael Luciuk", "mrl280@usask.ca",
                                  "University of Saskatchewan", format="hdf5")
            i = i + 1
