import numpy as np
from scipy.integrate import solve_ivp
def pulsazione(m,k):
  """
  Calcola la pulsazione naturale (frequenza angolare) di un oscillatore armonico.

  Parametri:
  m (float): Massa del sistema.
  k (float): Costante elastica della molla.

  Ritorna:
  float: La pulsazione naturale w0.
  """
  w0=np.sqrt(k/m)
  return w0
def rapporto_smorzamento(m,k,c):
  """
  Calcola il rapporto di smorzamento (zeta) dell'oscillatore.

  Parametri:
  m (float): Massa del sistema.
  k (float): Costante elastica della molla.
  c (float): Coefficiente di smorzamento viscoso.

  Ritorna:
  float: Il rapporto di smorzamento zeta.
  """
  zeta=c/(2*np.sqrt(k*m))
  return zeta
def regime(zeta):
  """
  Determina il tipo di regime di smorzamento in base al valore di zeta.

  Parametri:
  zeta (float): Il rapporto di smorzamento.

  Ritorna:
  str: Una stringa che identifica il regime ('sottosmorzato', 'critico' o 'sottosmorzato').
  """
  if np.isclose(zeta,1.0):
    regime_z="critico"
  elif zeta>1:
    regime_z="sovrasmorzato"
  else:
    regime_z="sottosmorzato"
  return regime_z
def energia(m,k,x,v):
  """
  Calcola l'energia totale (meccanica) del sistema in un dato istante.

  Parametri:
  m (float): Massa del sistema.
  k (float): Costante elastica della molla.
  x (float o array): Posizione/spostamento.
  v (float o array): Velocità.

  Ritorna:
  float o array: L'energia totale (somma di energia cinetica ed elastica).
  """
  E=0.5*m*v**2+0.5*k*x**2
  return E
##########aggiornamento es 2
def x_esatta(t,x0,m,k,c):
  """
  Calcola la soluzione analitica esatta dello spostamento x(t) per un oscillatore armonico smorzato
  con condizioni iniziali x(0) = x0 e v(0) = 0.

  Parametri:
  t (float o array): Istante o istanti di tempo in cui calcolare lo spostamento.
  x0 (float): Posizione iniziale all'istante t = 0.
  m (float): Massa del sistema.
  k (float): Costante elastica della molla.
  c (float): Coefficiente di smorzamento viscoso.

  Ritorna:
  float o array: Lo spostamento x(t) calcolato per i parametri e i tempi forniti.
  """
  w0=pulsazione(m,k)
  zeta=rapporto_smorzamento(m,k,c)
  if np.isclose(zeta,1.0):
    x=x0*np.exp(-w0*t)*(1+(w0*t))
  elif zeta <1:
    w_d=w0*(np.sqrt(1-zeta**2))
    x=x0*np.exp(-zeta*w0*t)*(np.cos(w_d*t)+(zeta*w0/w_d)*np.sin(w_d*t))
  else:
    w_star=w0*(np.sqrt(zeta**2-1))
    x=x0*np.exp(-zeta*w0*t)*(np.cosh(w_star*t)+(zeta*w0/w_star)*np.sinh(w_star*t))
  return x
def derivate(t,y,m,k,c):
  """
  Calcola le derivate temporali del vettore di stato [x, v] per l'equazione differenziale dell'oscillatore.

  Parametri:
  t (float): Istante di tempo corrente (richiesto da solve_ivp).
  y (list o array): Vettore di stato corrente [x, v], dove x è la posizione e v è la velocità.
  m (float): Massa del sistema.
  k (float): Costante elastica della molla.
  c (float): Coefficiente di smorzamento viscoso.

  Ritorna:
  list: Le derivate [dxdt, dvdt], ovvero la velocità e l'accelerazione istantanee.
  """
  x,v=y[0],y[1]
  dvdt=-(k/m)*x-(c/m)*v
  dxdt=v
  return [dxdt,dvdt]
def simula_scipy(m,k,x0,T_sim,c=0.0,v0=0.0,rtol=1e-9,atol=1e-12,n_punti=6001):
  """
  funzione che implementa il metodo di solve_ivp di scipy
  dati:
  massa(m)
  costante di elasticità(k)
  posizione iniziale(x0)
  durata(T_sim)
  coeff. di attrito viscoso(c)
  tolleranze (rtol,atol)
  numero di punti(n_punti)

  restituisce:
  soluzione problema iniziale con solve_ivp (y=[x,v])
  array tempo (t)

  """
  t_eval=np.linspace(0,T_sim,n_punti)
  sol=solve_ivp(derivate,(0,T_sim), [x0,v0], t_eval=t_eval, rtol=rtol, atol=atol, args=(m,k,c))
  t=t_eval
  if sol.success==False:
    raise RuntimeError(sol.message)
  return t,sol.y
def tempo_assestamento(t,x,soglia=0.001):
  """
  Trova l'indice dell'ultimo istante di tempo in cui lo spostamento dell'oscillatore supera o è uguale a una determinata soglia.

  Parametri:
  t (array): Vettore dei tempi.
  x (array): Vettore degli spostamenti nel tempo.
  soglia (float, opzionale): Valore di soglia assoluto per definire l'assestamento. Il valore di default è 0.001.

  Ritorna:
  int o float: L'ultimo istante di tempo in cui la soglia viene superata, oppure 0.0 se non viene mai superata.
  """
  superamenti = np.where(np.abs(x) >= soglia)[0]
  if len(superamenti) == 0:
    return 0.0
  t_assestamento = t[superamenti[-1]]
  return t_assestamento
