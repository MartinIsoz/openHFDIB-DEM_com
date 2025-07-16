#!/usr/bin/python

import time
start_time = time.time()
import re                                                       #regexp
import numpy as np
import os
from math import ceil
from math import floor
import sys

from scipy.integrate import odeint

import matplotlib.pyplot as plt

def surrogate(radius,g,rhoF,muF,rhoP,s0=0.0,v0=0.0,tE=1.0):
    k = 10e6                                                            #hardcoded, not used
    c = 100
    
    mp = (np.pi*radius**3*4/3)*rhoP
    mf = (np.pi*radius**3*4/3)*rhoF
    
    t = np.linspace(0,tE,1000)
    x0 = [v0, s0]
    
    resX,outMsg = odeint(model, x0, t, args=(mp,mf,radius,rhoF,g,k,c,muF,),full_output=True)
    vel = resX[:,0]
    s = resX[:,1]
    f_drag = comp_drag(vel,radius,mp,mf,g,[muF,rhoF])
    resX = np.hstack((resX,np.expand_dims(f_drag,axis=1)))
    
    return resX,t
    

def model(x,t,mp,mf,radius,rhoF,g,k,c,muF):
    Fhit = 0
    if x[1] <= 0:
        Fhit = k*x[1] + c*x[0]
    
    v = x[0]
    gravity = (mp-mf)*g
    pars= [muF,rhoF]
    drag = 0.5*np.pi*radius**2*rhoF*np.abs(v)*v*Cd(Re(v,radius,pars))
    
    # ~ print(t,v)
    
    dv_dx = (gravity-drag-Fhit)/(mp+0.0*mf)
    ds_dx = v
    return [dv_dx, ds_dx]
    
    
def Re(v,radius,pars):
    muF,rhoF = pars
    nu = muF/rhoF
    v = np.abs(v)
    v += 1e-12
    return 2*radius*v/nu
    
def comp_drag(v,radius,mp,mf,g,pars):
    muF,rhoF = pars
    gravity = (mp-mf)*g
    return (gravity - 0.5*np.pi*radius**2*rhoF*np.multiply(np.abs(v),v)*Cd(Re(v,radius,pars)))

def Cd(Re):
    return (24./Re)*(1.+0.27*Re)**0.43+0.47*(1.-np.exp(-0.04*Re**0.38)) # Cheng 2009 Comparison of formulas for drag...
    #return (24./Re)*(1.+0.15*Re**0.657)+0.407/(1.+8710*(1/Re))
    #return (777*((669806/875)+(114976/1155)*Re+(707/1380)*Re**2))/(646*Re*((32869/952)+(924/643)*Re+(1/385718)*Re**2))  # Empirical formula Mikhailov 2013 The drag coefficient of a sphere...

def readData(data,lookedForString,timeString,debuging=False):
    
    if(isinstance(lookedForString, str)):
        lookedForString = [lookedForString]
    elif(isinstance(lookedForString, list)):
        lookedForString = lookedForString
    
    if debuging:
        print("lookedForString: ", lookedForString)

    outerList = [[] for i in range(len(lookedForString)+1)]
    for i in range(len(data)):
        fInd = data[i].find(timeString)

        if fInd==0:
            outerList[-1].append([float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", data[i][fInd+len(timeString)::])])
            for pos in range(len(lookedForString)):
                string = lookedForString[pos]
                if len(outerList[-1]) > 0:
                    for k in range(1,2000):
                        if i-k > 0:
                            fInd = data[i - k].find(string)    
                            if fInd==0:
                                outerList[pos].append([float(strToVal) for strToVal in re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", data[i - k][fInd+len(string)::])])
                                break
                        else:
                            break        

    simTime = outerList[-1]
    time = [val[0] for val in simTime]
    dataLists = outerList[:-1]
    
    exportData = [[] for i in range(len(lookedForString))]

    if debuging:
        print("len(dataLists): ", len(dataLists))
        for pos in range(len(dataLists)):
            print("len(dataLists[",pos,"]): ", len(dataLists[pos]))
        print("len(simTime): ", len(simTime))


    for pos in range(len(dataLists)):
        dataList = dataLists[pos]
        if(len(dataList) > 0):
            for data in dataList:
                if(len(data) == 3):
                    exportData[pos].append(np.array(data))
                else:
                    exportData[pos].append(data)

    if debuging:
        print("len(exportData): ", len(exportData))
        for pos in range(len(exportData)):
            print("len(exportData[",pos,"]): ", len(exportData[pos]))
    
    return exportData,time

def writeOutputFile(data, headers, outputFile, collWidth=20,nFrom=0):
    # Define column width
    collumnWidth = collWidth

    with open(outputFile, "w") as f:
        # Print header
        f.write("".join(h.ljust(collumnWidth) for h in headers) + "\n")
        for row in data[nFrom::]:
            paddedRow = row + [""] * (len(headers) - len(row))
            f.write("".join(f"{str(val):<{collumnWidth}}" for val in paddedRow) + "\n")

def readDataFromFile(fileName,timeTreshold=1e6):
    with open(fileName, 'r') as file:
        data = file.readlines()
        nCol = len(data[0].split())
        # print(" in file: ", fileName, " nCol: ", nCol)
        readData = [[] for i in range(nCol)]
        for line in data[1:]:
            line = line.strip()
            for i in range(nCol):
                try:
                    # print(line.split()[i])
                    readData[i].append(float(line.split()[i]))
                except:
                    readData[i].append(line.split()[i])
            if(float(line.split()[0]) > timeTreshold):
                break
        
        file.close()
        return readData
        
def list_folders(directory):
    return [name for name in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, name))]
        
def create_folder(folder_name):
    # ~ os.makedirs(folder_name, exist_ok=True)

    if(not os.path.exists(folder_name)):
        os.makedirs(folder_name)
    else:
        os.system('rm -r ' + folder_name)
        os.makedirs(folder_name)
    
def read_dict_value(filename,var_nm):
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('//') or line.startswith('/*') or line.startswith('*') or line.startswith('//~'):
                continue
            if var_nm in line:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == var_nm:
                    return float(parts[-1].rstrip(';'))                 # the value is the last item, assume the line ends with ";"
    return None                                                         # if not found
    
def read_dict_vector(filename,var_nm):
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('//') or line.startswith('/*') or line.startswith('*') or line.startswith('//~'):
                continue
            if var_nm in line:
                parts = line.split("(")[-1]
                parts = parts.split(");")[0]
                parts = parts.split(" ")
                if len(parts) >= 1:
                    return np.array([float(val) for val in parts])
    return None 

def main():
    
    inDir  = "./20_testCases/"
    outDir = "./ZZ_post_processed_results/"
    inFile = ['log.pimpleHFDIBFoam']
    plotData = True
    saveFig = True
    
    testIdStr = "settling_"                                             #take into account only these
    skipIdStr = "_clean"                                                #skip these
    
    # -- parameters that are assumed constant
    g = -9.81

    soughtDict = {
        '-- body 0 linear velocity      : (': 'velocity',
        '-- body 0 CoM                  : (': 'position',
        '-- body 0 Acting Force    : (':'force',
    }

    soughtStrings = list(soughtDict.keys())
    
    compFolders = list_folders(inDir)
            
    for folder in compFolders:
        
        if folder.find(testIdStr) < 0: continue
        if folder.find(skipIdStr) >= 0: continue
        
        print("Working on the case: " + folder)

        baseFolder = outDir + folder + "/"
        
        create_folder(baseFolder)
    
        exportData = []
        simulationTime = []
    
        with open(inDir + folder + "/" + inFile[0], 'r') as file:
            data = file.readlines()
            exportedData, simulationTime = readData(data,soughtStrings,'Time = ')
            file.close()
        
        # -- load the simulation settings for the surrogate model to create comparison
        nuF = read_dict_value(inDir + folder + "/" + "constant/transportProperties", "nu")
        rhoF= read_dict_value(inDir + folder + "/" + "constant/transportProperties", "rho")
        muF = nuF*rhoF
        
        rhoP= read_dict_value(inDir + folder + "/" + "constant/HFDIBDEMDict", "rho")   #works only for single material
        radius= read_dict_value(inDir + folder + "/" + "constant/HFDIBDEMDict", "radius")#works only for single object
        
        startPos = read_dict_vector(inDir + folder + "/" + "constant/HFDIBDEMDict","startPosition")#works only for single object
        s0     = startPos[1]
        
        tEnd= read_dict_value(inDir + folder + "/" + "system/controlDict", "endTime")#works only for single object
            
        surrRes,surrT = surrogate(radius,g,rhoF,muF,rhoP,s0,0.0,tEnd)
            
        # exporting the data 
        print("len(exportedData)\t: ", len(exportedData))
        print("len(simulationTime)\t: ", len(simulationTime))
        for pos in range(len(soughtStrings)):
    
            dataSet = exportedData[pos]
            print("len("+soughtDict[soughtStrings[pos]]+")\t: ", len(dataSet))
    
            bodyId = re.findall(r'\d+',soughtStrings[pos])[0]
    
            maxLength = min(len(dataSet),len(simulationTime))
    
            Header = ['t']
            Rows = []
            try:
                if(len(dataSet[0]) == 3):
                    Header.append('x')
                    Header.append('y')
                    Header.append('z')
                else:   
                    Header.append('data')                
            except:
                Header.append('data')
            
            for i in range(0,maxLength):
                row = []
                row.append(simulationTime[i])
                try:
                    for data in dataSet[i]:
                        row.append(data)
                except:
                    row.append(dataSet[i])
                Rows.append(row)
            
    
            fileName = baseFolder+'Data_' + str(bodyId) +'_'+ soughtDict[soughtStrings[pos]]+ '.txt'
            writeOutputFile(Rows, Header, fileName, 20)
    
        # Plotting the data
        if(plotData):
            ref = len(list(soughtDict.keys()))
            b = ceil(ref/2)
            focusOn = 'y'
            fig, axes = plt.subplots(2, b, figsize=(16, 10))  # 3 rows, 2 columns
            axes = axes.flatten()
            for pos in range(len(soughtStrings)):
                fileName = baseFolder+'Data_' + str(bodyId) +'_'+ soughtDict[soughtStrings[pos]]+ '.txt'
                data =readDataFromFile(fileName,0.2)
                dataToPlot = []
                time = []
                if(len(data) ==2):
                    time = data[0]
                    dataToPlot = data[1]
                elif(len(data) ==4):
                    time = data[0]
                    if(focusOn == 'x'):
                        dataToPlot = data[1]
                    elif(focusOn == 'y'):
                        dataToPlot = data[2]
                    elif(focusOn == 'z'):
                        dataToPlot = data[3]
                    else:
                        print("focusOn is not defined")
                        dataToPlot = data[2]
    
                axes[pos].plot(time, dataToPlot,"bo",markersize=1)
                axes[pos].plot(surrT,surrRes[:,pos],'red')
                axes[pos].set_xlabel('Time')
                axes[pos].set_title(soughtDict[soughtStrings[pos]])
    
            plt.suptitle(folder + "; nu = %10.2e (m2/s); rhoS/rhoF = %5.3f"%(nuF,rhoP/rhoF))
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            if saveFig:
                plt.savefig(outDir + folder + ".png")
            else:
                plt.show()
    

if __name__ == "__main__":
    main()
