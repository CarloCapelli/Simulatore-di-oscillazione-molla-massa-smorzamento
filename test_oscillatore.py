import numpy as np
import oscillatore
from oscillatore import pulsazione,rapporto_smorzamento,regime,energia,x_esatta,derivate,simula_scipy,tempo_assestamento
m,k,x0=2,200,0.1
def test_pulsazione():
  w=pulsazione(m,k)
  assert np.isclose(w,10)
def test_rapporto():
  zeta=rapporto_smorzamento(m,k,4)
  assert np.isclose(zeta,0.1)
  zeta_2=rapporto_smorzamento(m,k,40)
  assert np.isclose(zeta_2,1)
def test_regime():
  regime_1=regime(0.1)
  assert regime_1=="sottosmorzato"
  regime_2=regime(1)
  assert regime_2=="critico"
  regime_3=regime(2)
  assert regime_3=="sovrasmorzato"
def test_x0():
  for c in [4,40,80]:
    x=x_esatta(0.0,x0,m,k,c)
    assert np.isclose(x,x0)
def test_sim_contro_esatta():
  for c in [4,40,80]:
    t,sol=simula_scipy(m,k,x0,6,c=c)
    x_sol=sol[0]
    x_esatt=x_esatta(t,x0,m,k,c)
    assert np.max(np.abs(x_sol-x_esatt))<1e-8
def test_energia():
  t,sol=simula_scipy(m,k,x0,6,4)
  x=sol[0]
  v=sol[1]
  E=energia(m,k,x,v)
  assert np.all(np.diff(E)<1e-9)
def test_regime_vicino_a_1():
  tst=regime(1+1e-10)
  assert tst=="critico"
def test_assestamento_mai_sopra_soglia():
  t=np.linspace(0,1,11)
  x=np.zeros(11)
  tst=tempo_assestamento(t,x)
  assert tst==0.0
def test_assestamento_istante():
  t=np.linspace(0,0.5,6)
  x=np.array([1,1,1,0,0,0])
  t_assestamento=tempo_assestamento(t,x)
  assert t_assestamento==0.2
def test_x_esatta_critico():
  c,t=40,0.5
  x=x_esatta(t,x0,m,k,c)
  assert np.isclose(x,0.0040428)
