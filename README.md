# Simulatore-di-oscillazione-molla-massa-smorzamento
Simulatore comportamento di oscillazione di una molla con smorzamento. Destinato a studenti di fisica e meccanica che desiderano capire e simulare i movimenti di oscillazione smorzati

#Installazione
Clona il repository sul tuo computer locale e installa tutte le dipendenze necessarie tramite il file `requirements.txt` con i seguenti comandi:

``bash
git clone (https://github.com/CarloCapelli/Simulatore-di-oscillazione-molla-massa-smorzamento.git)
cd Simulatore-di-oscillazione-molla-massa-smorzamento
pip install -r requirements.txt

# Requisiti

Per eseguire questo progetto è necessario aver installato **Python 3.8** o versioni successive.

# Esempio d'uso
import oscillatore

import importlib

importlib.reload(oscillatore)

from oscillatore import pulsazione,rapporto_smorzamento,regime,energia,x_esatta,derivate,simula_scipy,tempo_assestamento

m,k,x0=2,200,0.1

c=[4,40,80]

t_sim=np.linspace(0,6,6001)

w0=pulsazione(m,k)

print(f"w0={w0}")

for c_sim in c:

  zeta=rapporto_smorzamento(m,k,c_sim)
  
  regime_sim=regime(zeta)
  
  print(f"c={c_sim}, zeta={zeta}, regime={regime_sim}")
  
  x=x_esatta(t_sim,x0,m,k,c_sim)
  
  assert np.isclose(x[0],x0), "posizione iniziale non corretta"
  
  t,sol=simula_scipy(m,k,x0,6,c_sim)
  
  tempo_ax=tempo_assestamento(t,sol[0])
  
  errore_max=np.max(np.abs(sol[0]-x))
  
  print(f"tempo assestamento = {tempo_ax} s")
  
  print(f"errore massimo = {errore_max:.2e}")

output atteso:

w0=10.0

c=4, zeta=0.1, regime=sottosmorzato

tempo assestamento = 4.48 s

errore massimo = 4.96e-11

c=40, zeta=1.0, regime=critico

tempo assestamento = 0.663 s

errore massimo = 6.15e-12

c=80, zeta=2.0, regime=sovrasmorzato

tempo assestamento = 1.746 s

errore massimo = 1.16e-11

# Come eseguire i test

I test automatici sono implementati con `pytest`. Per eseguirli e verificare il corretto funzionamento dell'oscillatore, assicurati di essere nella cartella principale del progetto e lancia il seguente comando nel terminale:

``bash
pytest -v

# Limiti

Questo progetto implementa una simulazione ideale che presenta alcune limitazioni fisiche e numeriche da considerare:

* **Parametri solo positivi:** Il modello fisico assume che massa ($m > 0$), costante elastica della molla ($k > 0$) e coefficiente di smorzamento ($c \ge 0$) siano strettamente positivi o non negativi. Valori negativi porterebbero a comportamenti fisicamente instabili.
* **Smorzamento viscoso lineare:** Si assume che la forza di smorzamento sia direttamente proporzionale alla velocità ($F_v = -c \cdot v$). Non sono considerati attriti di tipo coulombiano (attrito secco/radente) o quadratico.
* **Comportamento elastico lineare:** La molla segue rigorosamente la legge di Hooke lineare ($F_e = -k \cdot x$). Non vengono gestite deformazioni plastiche o molle non lineari.
* **Passo d'integrazione:** Per simulazioni estremamente rapide o con costanti di smorzamento molto elevate (sistemi fortemente rigidi), potrebbe essere necessario regolare manualmente le tolleranze di integrazione (`rtol` e `atol`) o aumentare il numero di punti temporali per evitare errori di troncamento numerico.
