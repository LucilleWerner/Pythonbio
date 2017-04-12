#programma voor berekenen van het percentage cysteine, tryptofaan, de meest voorkomende aminozuren,
#de minst voorkomende aminozuren, de hydrofiele aminozuren en hydrophobe aminozuren.
#LucilleWerner
#08-03-16

from collections import Counter

def main():
    a = open("HIV1.txt")
    b = open("HIV2.txt")
    c = open("SIVmnd2.txt")
    files = [a,b,c]

    names = ["HIV1", "HIV2", "SIVmnd2"]
    genen = ["externe proteinen", "interne proteinen"]

    #per virus worden de envelop proteinen (externe proteinen) gescheiden van de rest (interne proteinen)
    #de namen van de virussen worden ook meegenomen in deze forloop, ze worden later synchroon met de resultaten van de virussen geprint
    
    for bestand, name in zip(files, names):
        proteins = []
        env, rest = fastaconvert(bestand)
        proteins.append(env)
        proteins.append(rest)
        #print("env percentage ",len(env)/(len(env)+len(rest)))
        #print("rest percentage ",len(rest)/(len(env)+len(rest)))
        #per sets proteins (envelop en rest) worden alle berekeningen van het script uitgevoerd
        #het resultaat is dat alle berekeningen per virus per proteinen (envelop en rest) worden uitgevoerd en geretourneerd
        
        for protein, gen in zip(proteins, genen):
            cysper, tryper = cysteine_tryptofaan(protein)
            maxfreqs, minfreqs = minmaxfreqs(protein) 
            phoper, phiper, inper = hydrophile_hydrophobe(protein)
            
            print("Het cysteinepercentage van", gen, "van", name, "is", ("%.1f" % cysper), "%" )
            print("Het tryptofaanpercentage van", gen, "van", name, "is", ("%.1f" % tryper), "%" )
            print("De meestvoorkomende aminozuren van", gen, "van", name, "zijn", maxfreqs[0][0],":",maxfreqs[0][1],"%", maxfreqs[1][0],":", maxfreqs[1][1],"%", maxfreqs[2][0],":", maxfreqs[2][1],"%")
            print("De minstvoorkomende aminozuren van", gen, "van",  name, "zijn",minfreqs[0][0],":",minfreqs[0][1],"%", minfreqs[1][0],":", minfreqs[1][1],"%", minfreqs[2][0],":", minfreqs[2][1],"%")
            print("Het percentage hydrohobe aminozuren van", gen, "van", name, "is", ("%.2f" % phoper))
            print("Het percentage hydrofiele aminozuren van", gen, "van", name, "is", ("%.2f" % phiper), "%" )
            print("Het percentage aminozuren tussenin van", gen, "van", name, "is", ("%.2f" % inper), " %" )

            
def fastaconvert(bestand):
    bestand = bestand.readlines()
    
    #wanneer ">" in de regel zit worden de regels toegevoegd aan rest = []
    #wanneer "env" en niet ">" in de regel zit, wordt deze aan env = [] toegevoegd
    #wanneer hierna ">" gevonden wordt, is reread True en envread False, hierdoor worden de overige regels aan rest toegevoegd
    #alvorens het toevoegen aan lijsten worden de newlines eruit gefilterd en grote letters worden kleine letters
    #uiteindelijk worden er strings geretourneerd door "".join()
    
    envread = False
    reread = False
    env = []
    rest = []
    for regel in bestand:
        if ">" in regel:
            reread  = True
            envread  = False
        if 'env' in regel:
            envread = True
            reread = False
        if envread and ">" not in regel:
            regel = regel.replace("\n","")
            regel = regel.lower()
            env.append(regel)
        if reread and ">" not in regel:
            regel = regel.replace("\n","")
            regel = regel.lower()
            rest.append(regel)

    env = "".join(env)        
    rest = "".join(rest)
    
    return env, rest


def cysteine_tryptofaan(protein):

    #cysteine en tryptofaan percentages worden per protein (envelop en rest) berekend
    #unit is de lengte van het betreffende proteine
    
    unit = len(protein) 
    cyscount = protein.count("c")
    trypcount = protein.count("w")
    cysper = (cyscount/unit)*100
    tryper = (trypcount/unit)*100
    
    return cysper, tryper
    

def minmaxfreqs(protein):

    #de aminozuren met de maximale frequentie worden met de functie Counter berekend, hierover wordt .most_common() gecast
    #de aminozuren met de minimale frequentie worden op dezelfde manier berekend, alleen dan met slicing zodat de minst voorkomende berekend worden
    #in de variablen maxfreq en minfreq zit een dictionary die geconverteerd is naar een lijst vb: [('k', 183) ('l', 183)]
    #deze meest en minst voorkomende aminozuren worden samen met hun naam (1 letter) en frequentie toegevoegd aan aan maxlist en minlist
    
    total = len(protein) 
    maxfreq = list(Counter(protein).most_common(3))
    minfreq = Counter(protein).most_common(3)[-total:] 
    maxlist = []
    for i in range(3):
        freq = round((maxfreq[i][1]/total)*100, 2)
        maxlist.append([maxfreq[i][0], freq])
    minlist = []
    for i in range(3):
        freq = round((minfreq[i][1]/total)*100, 2)
        minlist.append([minfreq[i][0], freq])
        
    return maxlist, minlist
        
    
def hydrophile_hydrophobe(protein):

    #alle aminozuren die hydrofoob, hydrofiel en tussenin zijn staan in verschillende lijsten
    #deze lijsten staan alle drie genest in de lijst "properties", over deze lijst wordt geloopt
    #vervolgens worden over alle aminozuren in deze lijsten (phoob, phile en inbetween) geloopt
    #deze resultaten worden aan een nieuwe lijst toegevoegd (prop)
    #voor alle eigenschappen wordt een percentage berekend: delen door de totale lengte van protein en vermenigvuldigen met 100
    
    totaal = len(protein) 
    phoob = ["l", "i", "f", "v", "m", "w", "c", "a", "p"]
    phile = ["r", "d", "k", "e", "h", "n", "q"]
    inbetween = ["g", "s", "t", "y"]
    properties = [phoob, phile, inbetween]
    prop = []
    for eigenschap in properties:
        count = 0
        for amino in eigenschap:
            count += protein.count(amino)
        prop.append(count)
    
    phoper = prop[0] 
    phiper = prop[1]
    inper = prop[2]
    phoper = (phoper/totaal)*100
    phiper = (phiper/totaal)*100
    inper = (inper/totaal)*100
    return phoper, phiper, inper


main()
