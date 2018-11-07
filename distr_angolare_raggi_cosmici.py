"""
Giulio Grechi - Matteo Barbetti
Laboratorio di fisica subnucleare, esperienza rivelatore Adamo

distribuzione angolare raggi cosmici rivelati
"""

"""
parametri per l'allineamento e rototraslazione del piano centrale
    EXT PARAMETER                                   STEP         FIRST   
    NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
     1     dx         104.41       0.16200       0.18426E-05   -77.312    
     2     dy        -378.16       0.23979       0.46097E-05    7.2879    
     3     dz        -337.68       0.77977       0.59722E-05   -27.504    
     4     phi        12.848        15.837       0.17420E-04   -8.0738    
     5     psi       -310.35        4.6241       0.50209E-05    29.750    
     6     tet       -2416.7        20.414       0.19572E-04    9.8478   
"""
import copy
import time
import ROOT
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import numpy as np
from scipy.stats.stats import pearsonr

"""
scelta visualizzazione in gradi "deg" o in radianti "rad"
scelta del binning dei due istogrammi per theta e phi
scelta degli eventi di cui fare il grafico tridimensionale (intervallo)
"""
rad_deg = "rad" #in deg non è implementata correttamente l'accettanza
binning_n_theta = 80
binning_n_phi = 80
n_ev_graph_min, n_ev_graph_max = 0,100
flag_his_2d = "COLS9"
nome_istogramma = "dist_ang_mu.root"

"""
"""
tempo_misura = 30*60*60
superficie = float(7.0*5.33)

"""
gli offset devono essere inseriti con il segno fornito dal programma di 
minimizzazione
"""
offset_coordinates = [104.41, -378.16, -337.68]
offset_angles = [12.848,-310.35,-2416.7]
z_coord = [0., -0.8e4, -1.6e4] #um

"""
in presenza di piu' file, vanno prima concatenati in uno solo
"""

nome_file = "track700n_totali.dat"

"""
opzioni possibili da passare al programma da riga di comando
"""
if "debug" in sys.argv: debug = True
else: debug = False

if "graph" in sys.argv: graph = True
else: graph = False

if "hist" in sys.argv: histogram = True
else: histogram = False

if "save" in sys.argv: save = True
else: save = False

if "rad" in sys.argv: rad_deg = "rad"
if "deg" in sys.argv: rad_deg = "deg" 

if "test" in sys.argv: test = True
else: test = False

srt_time = time.time()
random.seed()

"""
crea una lista leggendo un file, in cui ogni elemento e' a sua volta una lista 
contenente i valori numerici contenuti in una data riga
N         x1        x2        x3        y1        y2        y3
8.000000  27492.41  26772.95  26327.72  42574.12  45471.26  47433.26    
23.00000  29719.47  27803.59  26177.88  12598.04  7584.102  2111.649    
...
diventa
[
[8.0, 27492.41, 26772.95, 26327.72, 42574.12, 45471.26, 47433.26]
[23.0, 29719.47, 27803.59, 26177.88, 12598.04, 7584.102, 2111.649]
...]

"""
def readFileMeasure(file,first_line):
    datas = [
            [float(i) for i in x.strip().split()]
             for x in open(file,'r').readlines()[first_line:]
           ]
    if debug == True: print "[%2.4f]numero di eventi: %s"%(
                                          time.time() - srt_time, len(datas))
    return datas


"""
prende una lista di liste nella forma 
N     x1        x2        x3        y1        y2        y3
[[8.0, 27492.41, 26772.95, 26327.72, 42574.12, 45471.26, 47433.26]...]
e la trasforma in una del tipo
[ [ [x1 y1 z1],[x2 y2 z2],[x3 y3 z3]]...]
"""
def unpack(list_of_float):
    z_1 = z_coord[0]
    z_2 = z_coord[1]
    z_3 = z_coord[2]
    events_points = []
    for i in range(len(list_of_float)):
        events_points.append([
                              [list_of_float[i][1], list_of_float[i][4], z_1],
                              [list_of_float[i][2], list_of_float[i][5], z_2],
                              [list_of_float[i][3], list_of_float[i][6], z_3]
                              ])
    if debug == True: print "[%2.4f]creata lista di coordinate"%(
                                                        time.time() - srt_time)
    return events_points

"""
crea una nuova lista copiando quella in ingresso, del tipo 
[ [ [x1 y1 z1],[x2 y2 z2],[x3 y3 z3]]...]
le cui coordinate 2 sono spostate delle quantita' definita in offset
NB LE QUANTITa' IN OFFSET SONO SOMMATE!
"""
def shift(list_of_coordinates, offset):
    #nb: la lista non viene copiata ma ne viene cambiato solo il nome
    shifted_coordinates = list_of_coordinates 
    if debug == True: print "[%2.4f]evento non shiftato: \n%s"%(
                          time.time() - srt_time, shifted_coordinates[0][1])
    for i in range(len(shifted_coordinates)):
        for j in range(len(offset_coordinates)):
            shifted_coordinates[i][1][j] = (shifted_coordinates[i][1][j]
                                            + offset[j])
    if debug == True: print "[%2.4f]evento shiftato: \n%s"%(
                          time.time() - srt_time, shifted_coordinates[0][1])
    if debug == True: print "[%2.4f]applicato offset allineamento coordinate"%(
                                                       time.time() - srt_time,)
    return shifted_coordinates

"""
crea una nuova lista copiando quella in ingresso, del tipo 
[ [ [x1 y1 z1],[x2 y2 z2],[x3 y3 z3]]...]
le cui coordinate 2 sono ruotate
DA DEFINIRE
"""
def rotate(list_of_coordinates, angles):
    rotated_coordinates = list_of_coordinates
    if debug == True: 
        print "[%2.4f]se implementata, effettuata rotazione piano centrale"%(
                                                       time.time() - srt_time,)
    return rotated_coordinates

"""
accetta liste contenenti punti formattati dalla funzione unpack.
viene effettuato un fit per ogni tripletta di punti, sulla proiezione della 
retta nei due piani xy e xz, dai cui coefficienti angolari si possono ricavare
gli angoli theta e phi
ritorna in uscita due liste di uguale lunghezza contenenti i valori theta e phi

"""
def fit_linear(points_vector, angle):
    if debug == True: print "[%2.4f]fit delle tracce"%(time.time() - srt_time,)
    num_eventi = len(points_vector)
    result_phi = [0.0]*num_eventi
    result_theta = [0.0]*num_eventi
    fitFunction = ROOT.TF1("retta3d","[0]+[1]*x",0,1e5)
    for event_n in range(num_eventi):
        xy_tgraph = ROOT.TGraph()
        for index in range(3):
            xy_tgraph.SetPoint(index,
                            points_vector[event_n][index][0],
                            points_vector[event_n][index][1],)
        fit_xy = xy_tgraph.Fit(fitFunction,"SQ") #phi
        tan_result_xy = fit_xy.Parameters()[1]
        result_xy = np.arctan(tan_result_xy)
        result_phi[event_n] = result_xy
        control = (points_vector[event_n][2][1] 
                  - points_vector[event_n][0][1])
        if tan_result_xy < 0:
            result_phi[event_n] = result_phi[event_n] + 2*np.pi
        if control > 0 and tan_result_xy > 0 :
            result_phi[event_n] = result_phi[event_n] + np.pi
        if control < 0 and tan_result_xy < 0 :
            result_phi[event_n] = result_phi[event_n] - np.pi

        xz_tgraph = ROOT.TGraph()
        for index in range(3):
            xz_tgraph.SetPoint(index,
                            points_vector[event_n][index][0],
                            points_vector[event_n][index][2],)
        fit_xz = xz_tgraph.Fit(fitFunction,"SQ") #theta
        result_theta[event_n] = np.abs(np.arctan(
                                    1/(np.cos(result_phi[event_n])
                                    *fit_xz.Parameters()[1])))
    if test == True:
        for indx in range(len(result_phi)):
            rnd = random.random()*2*np.pi - np.pi
            result_phi[indx] = rnd
        if debug == True: 
            print "[%2.4f]generata distrubuzione uniforme casuale di phi"%(
                                                       time.time() - srt_time,)
        for indx in range(len(result_theta)):
            rnd = random.random()*np.pi*0.5
            result_theta[indx] = rnd
        if debug == True: 
            print "[%2.4f]generata distrubuzione uniforme casuale di theta"%(
                                                       time.time() - srt_time,)
        
    if debug == True:
        print "[%2.4f]fit effettuato"%(time.time() - srt_time,)
        corr = pearsonr(result_theta,result_phi)
        print "[%2.4f]correlazione theta phi = %1.4f, p-value = %0.2F"%(
                                     time.time() - srt_time, corr[0], corr[1])    
    if angle == "deg":
        return np.rad2deg(result_theta), np.rad2deg(result_phi)
    elif angle == "rad":
        return result_theta, result_phi
    else:
        print "errore, gradi o radianti?"
        sys.exit(0)

def distribution(vector_data_points, angle, n_bin_theta, n_bin_phi,flag):
    theta, phi, = fit_linear(vector_data_points, angle)
    
    if angle == "deg":
        min_theta, max_theta = 0, 90
        min_phi, max_phi = 0, 360
    else:
        min_theta, max_theta = 0, np.deg2rad(90)
        min_phi, max_phi = 0., np.deg2rad(360)
    if debug == True:
        print "[%2.4f]creazione e riempimento istogramma"%(
                                                       time.time() - srt_time,)
    histogram_angular_2d = ROOT.TH2F("Distribuzione angolare",
                                     "Distribuzione angolare RC",
                                     n_bin_theta,min_theta,max_theta,
                                     n_bin_phi,min_phi,max_phi,
                                     )

    fluxXaccettanza_hist_1d = ROOT.TH1F("flusso X accettanza",
                        "flusso X accettanza",
                        n_bin_theta,min_theta,max_theta,)
                        
    accettanza_hist_1d = ROOT.TH1F("accettanza",
                        "accettanza",
                        n_bin_theta,min_theta,max_theta,)
    for i in range(len(theta)):
        histogram_angular_2d.Fill(theta[i],phi[i])
        fluxXaccettanza_hist_1d.Fill(theta[i])
        accettanza_hist_1d.Fill(theta[i])

    function_0 = ROOT.TF1("sin(x)*cos(x)","sin(x)*cos(x)", 0, np.pi/2)
    factor = superficie * tempo_misura * 2 * np.pi
    fluxXaccettanza_hist_1d.Divide(function_0,factor)
    if debug == True:
        int_ottavo_ang_sol = fluxXaccettanza_hist_1d.Integral(0,29)
        print"[%2.4f]divisa distribuzione per i fattori geometrici e cost."%(
                                                      time.time() - srt_time,)
        print "integrale di un ottavo dell'angolo solido: %s"%int_ottavo_ang_sol
        print "proiezione su tutto l'angolo solido: %s m^-2 s^-1"%(int_ottavo_ang_sol*80000)
        
    canvas = ROOT.TCanvas()
    canvas.Divide(2,2)
    canvas.cd(1)
    histogram_angular_2d.GetXaxis().SetTitle("theta")
    histogram_angular_2d.GetYaxis().SetTitle("phi")
    histogram_angular_2d.Draw(flag)
    canvas.cd(2)
    histogram_angular_1d_theta = histogram_angular_2d.ProjectionX(
                                                        "Distribuzione Theta")
    histogram_angular_1d_theta.GetXaxis().SetTitle("theta")
    histogram_angular_1d_theta.GetYaxis().SetTitle("dN/d(theta)")
    histogram_angular_1d_theta.Draw()
    canvas.cd(3)
    histogram_angular_1d_phi = histogram_angular_2d.ProjectionY(
                                                          "Distribuzione Phi")
    histogram_angular_1d_phi.GetXaxis().SetTitle("phi")
    histogram_angular_1d_phi.GetYaxis().SetTitle("dN/d(phi)")    
    histogram_angular_1d_phi.Draw()
    canvas.cd(4)
    fluxXaccettanza_hist_1d.GetXaxis().SetTitle("theta")
    fluxXaccettanza_hist_1d.GetYaxis().SetTitle("1/(cm^2 s)")    
    #flux = ROOT.TF1("cos(x)^2","[0]*cos(x)^2", 0, np.pi/2)
    #fluxXaccettanza_hist_1d.Fit(flux,"","",0, 0.5)
    #flux.Draw()
    fluxXaccettanza_hist_1d.Draw()
    """
    canvas.cd(6)
    function_1 = ROOT.TF1("sin(x)*cos(x)^3","sin(x)*cos(x)*cos(x)*cos(x)", 0, np.pi/2)
    factor = superficie * tempo_misura * 2 * np.pi * 4 / np.pi
    accettanza_hist_1d.Divide(function_1, factor)
    accettanza_hist_1d.GetXaxis().SetTitle("theta")
    accettanza_hist_1d.GetYaxis().SetTitle("a.u.")    
    #fit_func = ROOT.TF1("cos(x)^p0","[1]*cos(x)^[0]", 0, np.pi/2)
    #accettanza_hist_1d.Fit(fit_func,"","",0, 1)
    accettanza_hist_1d.Draw()
    """
    if debug == True: print "[%2.4f]istogrammi costruiti"%(
                                                       time.time() - srt_time,)
    null = raw_input("premere per continuare")

    return histogram_angular_2d
"""
crea un grafico 3d disegnando le tracce di tutti gli eventi compresi 
tra min e max. richiede un invio per uscire dalla visualizzazione (raw_input)
"""
def graph_track(points_vec,min,max):
    if debug == True: print "[%2.4f]crezione del grafico 3d"%(
                                                       time.time() - srt_time,)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')
    ax.set_zlabel('Z [cm]')
    for ev_n in range(min,max,1):
        ax.plot([points_vec[ev_n][0][0]*1e-4,
                 points_vec[ev_n][1][0]*1e-4,
                 points_vec[ev_n][2][0]*1e-4],
                [points_vec[ev_n][0][1]*1e-4,
                 points_vec[ev_n][1][1]*1e-4,
                 points_vec[ev_n][2][1]*1e-4],
                [points_vec[ev_n][0][2]*1e-4,
                 points_vec[ev_n][1][2]*1e-4,
                 points_vec[ev_n][2][2]*1e-4],)
    if debug == True: print "[%2.4f]grafico creato"%(time.time() - srt_time,)
    plt.show()
    null = raw_input("premere per continuare")
    return True

"""
apre e legge il file "nome_file" a partire dal numero di riga inserito
creazione della lista con le coordinate e modifica della stessa lista
applicando offset e rotazione (nb: alla fine tutte le liste sono uguali, 
viene evitata la duplicazione per risparmiare tempo e memoria)
"""
def main():
    data_points = unpack(readFileMeasure(nome_file,0))
    data_points_shifted = shift(data_points, offset_coordinates)
    data_points_rotated = rotate(data_points_shifted, offset_angles)
    
    
    if graph == True:
        graph_track(data_points_shifted, 
                    n_ev_graph_min, n_ev_graph_max)
    if histogram == True:
        distribution(data_points_rotated, 
                     rad_deg, 
                     binning_n_theta, binning_n_phi,
                     flag=flag_his_2d,)
    
    return True
        

if __name__ == '__main__':
    main()
