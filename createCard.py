import ROOT
import sys
import os
import subprocess

#example: python createCard.py results_PHYS14
#example: python createCard.py results_PHYS14 hyp_highpt_excl_sr cards/card.txt

#then get expected limits with: combine -M Asymptotic results_PHYS14/card.txt --run expected --noFitAsimov

#to add more nuisances edit Process, writeOneCardFromProcesses and then set values in writeOneCard

class Process:
    def __init__(self, mycount, myname, myrootf, myplot):
        self.count = mycount
        self.name = myname
        self.rootf = myrootf
        self.plot = myplot
        self.fakes = "-"
        self.TTV = "-"
        self.WZ = "-"
    def rate(self): 
        f = ROOT.TFile(dir+"/"+self.rootf)
        if f.Get(self.plot): return f.Get(self.plot).Integral()
        else: 
            print self.plot+" not found in "+self.rootf
            return 0

#write card regardless of number of processes (but make sure signal is first in list)
def writeOneCardFromProcesses(dir, plot, output, processes):
    line = "---------------------------------------------------------------"
    bin = "SS"
    card = open(str(dir)+'/'+str(output), 'w')
    card.write("imax 1  number of channels \n")
    card.write("jmax *  number of backgrounds \n")
    card.write("kmax *  number of nuisance parameters \n")
    card.write(line+"\n")
    for process in processes:
        card.write("shapes "+process.name+" * "+dir+"/"+process.rootf+" "+plot+" "+plot+"\n")
    card.write("shapes data_obs * "+dir+"/TTWJets_histos.root "+plot+" "+plot+"\n")#dummy for now, please use --noFitAsimov option
    card.write(line+"\n")
    card.write("bin "+str(bin)+"\n")
    card.write("observation "+str(processes[1].rate())+"\n")#dummy for now, please use --noFitAsimov option
    card.write(line+"\n")
    #bin
    card.write("%-20s %-5s " % ("bin",""))
    for process in processes: card.write("%-10s " % (bin))
    card.write("\n")
    #process count
    card.write("%-20s %-5s " % ("process",""))
    for process in processes: card.write("%-10s " % (process.count))
    card.write("\n")
    #process name
    card.write("%-20s %-5s " % ("process",""))
    for process in processes: card.write("%-10s " % (process.name))
    card.write("\n")
    #process rate
    card.write("%-20s %-5s " % ("rate",""))
    for process in processes: card.write("%-10.3f " % (process.rate()))
    card.write("\n")
    #separate
    card.write(line+"\n")
    #nuisance fakes
    card.write("%-20s %-5s " % ("fakes","lnN"))
    for process in processes: card.write("%-10s " % (process.fakes))
    card.write("\n")
    #nuisance TTV
    card.write("%-20s %-5s " % ("TTV","lnN"))
    for process in processes: card.write("%-10s " % (process.TTV))
    card.write("\n")
    #nuisance WZ
    card.write("%-20s %-5s " % ("WZ","lnN"))
    for process in processes: card.write("%-10s " % (process.WZ))
    card.write("\n")
    return

def writeOneCard(dir, plot, output):
    #define processes (signal first)
    T1tttt = Process(0,"T1tttt","T1ttttG1200_histos.root",plot)
    TTW = Process(1,"TTW","TTWJets_histos.root",plot)
    TTZ = Process(2,"TTZ","TTZJets_histos.root",plot)
    WZ  = Process(3,"WZ","WZJets_histos.root",plot)
    ttbar = Process(4,"ttbar","ttbar_histos.root",plot)
    #overwrite nuisances
    TTW.TTV = "1.2"
    TTZ.TTV = "1.2"
    WZ.WZ = "1.2"
    ttbar.fakes = "1.5"
    #fill list of processes    
    processes = []
    processes.append(T1tttt)
    processes.append(TTW)
    processes.append(TTZ)
    processes.append(WZ)
    processes.append(ttbar)
    #create it
    writeOneCardFromProcesses(dir, plot, output, processes )
    return

def writeAllCards(dir):
    writeOneCard(dir, "hyp_highpt_sr", "card-hi-hi.txt" )
    writeOneCard(dir, "hyp_lowpt_sr", "card-hi-low.txt" )
    writeOneCard(dir, "hyp_verylowpt_sr", "card-low-low.txt" )
    olddir = os.getcwd()
    os.chdir(dir)
    #os.system("python combineCards.py card-hi-hi.txt card-hi-low.txt card-low-low.txt >& card-all.txt")
    f = open('card-all.txt', 'wb')
    subprocess.call(["combineCards.py","card-hi-hi.txt","card-hi-low.txt","card-low-low.txt"],stdout=f)
    os.chdir(olddir)

#main body
if len(sys.argv)==2:
    dir = sys.argv[1]
    writeAllCards( dir )
elif len(sys.argv)==4:
    dir = sys.argv[1]
    plot = sys.argv[2]
    output = sys.argv[3]
    writeOneCard( dir, plot, output )
else: print "number of arguments not supported"
