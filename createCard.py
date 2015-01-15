import ROOT

dir = "results_PHYS14_new_noptrel"

class Process:
    def __init__(self, mycount, myname, myrootf):
        self.count = mycount
        self.name = myname
        self.rootf = myrootf
        self.fakes = "-"
        self.TTV = "-"
        self.WZ = "-"
    def rate(self): 
        f = ROOT.TFile(dir+"/"+self.rootf)
        return f.Get("hyp_highpthtmetmt_sr").Integral()

T1ttttG1500 = Process(0,"T1tttt","T1ttttG1500_histos.root")
TTW = Process(1,"TTW","TTWJets_histos.root")
TTZ = Process(2,"TTZ","TTZJets_histos.root")
WZ  = Process(3,"WZ","WZJets_histos.root")
ttbar = Process(4,"ttbar","ttbar_histos.root")

TTW.TTV = "1.2"
TTZ.TTV = "1.2"
WZ.WZ = "1.2"
ttbar.fakes = "1.5"

processes = []
processes.append(T1ttttG1500)
processes.append(TTW)
processes.append(TTZ)
processes.append(WZ)
processes.append(ttbar)

line = "---------------------------------------------------------------"
bin = "SS"

card = open('card.txt', 'w')

card.write("imax 1  number of channels \n")
card.write("jmax *  number of backgrounds \n")
card.write("kmax *  number of nuisance parameters \n")
card.write(line+"\n")
for process in processes:
    card.write("shapes "+process.name+" * "+dir+"/"+process.rootf+" hyp_highpthtmetmt_sr"+" hyp_highpthtmetmt_sr"+"\n")
card.write(line+"\n")
card.write("bin "+str(bin)+"\n")
card.write("observation "+str(0)+"\n")
card.write(line+"\n")
card.write("%-20s %-5s %-10s %-10s %-10s %-10s %-10s \n" % ("bin","",bin,bin,bin,bin,bin))
card.write("%-20s %-5s %-10s %-10s %-10s %-10s %-10s \n" % ("process","",T1ttttG1500.count,TTW.count,TTZ.count,WZ.count,ttbar.count))
card.write("%-20s %-5s %-10s %-10s %-10s %-10s %-10s \n" % ("process","",T1ttttG1500.name,TTW.name,TTZ.name,WZ.name,ttbar.name))
card.write("%-20s %-5s %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f \n" % ("rate","",T1ttttG1500.rate(),TTW.rate(),TTZ.rate(),WZ.rate(),ttbar.rate()))
card.write(line+"\n")
card.write("%-20s %-5s %-10s %-10s %-10s %-10s %-10s \n" % ("fakes","lnN",T1ttttG1500.fakes,TTW.fakes,TTZ.fakes,WZ.fakes,ttbar.fakes))
card.write("%-20s %-5s %-10s %-10s %-10s %-10s %-10s \n" % ("TTV","lnN",T1ttttG1500.TTV,TTW.TTV,TTZ.TTV,WZ.TTV,ttbar.TTV))
card.write("%-20s %-5s %-10s %-10s %-10s %-10s %-10s \n" % ("WZ","lnN",T1ttttG1500.WZ,TTW.WZ,TTZ.WZ,WZ.WZ,ttbar.WZ))

