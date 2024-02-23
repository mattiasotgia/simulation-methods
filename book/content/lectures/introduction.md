# Introduzione

## Metodi di generazione Monte Carlo (MC)

È un metodo utilizzato in un sacco di branche della fisica, dalla fisica delle particelle alla fisica della materia, è uno strumento che nasce negli anni '40, ed ha una continua evoluzione. 

:::{margin}
Le origini del Metodo di Monte Carlo risalgono alla metà degli anni 40 nell’ambito del Progetto Manhattan.
I formalizzatori del metodo sono Enrico Fermi, John von Neumann e Stanislaw Marcin Ulam, il nome Monte
Carlo fu inventato in seguito da Nicholas Constantine Metropolis in riferimento alla nota tradizione nei giochi
d’azzardo dello stato omonimo nel sud della Francia (e in particolare alla roulette in cui i numeri sono estratti a
caso).

Pare che l’idea venne a Ulam mentre era convalescente da una malattia. Per passare il tempo giocava a un
solitario con le carte e si chiese quale fosse la probabilità di completarlo correttamente. Le regole del solitario
rendevano difficile il calcolo di questa probabilità. Tuttavia pensò che si sarebbero potute simulare con un
computer molte disposizioni casuali di un mazzo di carte e controllare in quanti di questi casi si riusciva a
completare il solitario. In questo modo si poteva calcolare la probabilità di vincere al solitario in modo empirico.
Ulam ne parlò con Von Neumann il quale comprese il potenziale di questa intuizione e svilupp`o il primo
algoritmo che permetteva di generare numeri casuali con un computer.
:::

In fisica della materia, sopratutto in fisica teorica della materia, si parla di 
 - Calcolare proprietà di equilibrio di sistemi a molte particelle (N-body). In particolare eccetto il caso di interazioni armoniche, è praticamente impossibile ottenere una soluzione analitica, e quindi è obbligatorio utilizzare strumenti numerici e metodi numerici.
 - Processi stocastici (rndm walk, crescita di cristalli).

```{figure} ../../images/ATLAS-ttH-eventdisplay.png
---
figclass: margin
---
ATLAS Event display (for a $\mathrm{t\bar t H}$ event)
```

In fisica delle particelle (nuclear and sub-nuclear physics)
 - Simulation of the initial process (of which the out state is known by experimental data)
 - Simulazione dell'interazione con un rivelatore
 - Produzione di eventi identici a quelli veri, per creare maggiore statistica, e/o per progettare/ottiimizare/controllare il rivelatore.

## Numeri casuali

Centrale negli algoritmi di MC è la necessità di poter **generare** un numero casuale. In generale questo generatore deve avere alcune proprietà, alcune convensionali altre necessarie
 1. Deve generare in $[0, 1)$.
 2. Deve generare valori che non rispettano nessun pattern apparente ($x_n$ è scorrelato a $x_{n-1}$).
 3. Deve avere periodo di ripetizione estremamente lungo
 5. Deve essere sufficentemente veloce (sopratutto per utilizzo in particle physics, LHC)

Assumiamo adesso di possedere un generatore di numeri casuali, con le proprietà precedenti, e che sia tale per cui 

$$
    \eta \in [0, 1), \quad p(\eta) = 1, \quad \int_0^1 p(\eta) d \eta = 1.
$$

Dalla generazione di numeri casuali, ci interessa allora poter generare distribuzioni continue, e già abbiamo visto alcuni metodi
 - Inversione
    - Molto utilizzato per funzioni integramìbili,
    - Garantisce l'efficenza massima
 - Reiezione (senza campionamento, anche detto reiezione semplice)
    - Valido per qualsiasi funzione
    - Efficenza che però dipende dalla funzione che si sta considerando (quindi non è prevedibile a priori, ma è _model dependent_)
 - Reiezione (con campionamento/importance sampling)
    - Valido anche questo per tutte le funzioni
    - Rispetto all'algoritmo di reiezione può avere efficenza migliore, ma comunque resta _omdel dependent_

Esiste ancora un modo di generare, detto campionamento semplice.

### Campionamento semplice

È un metodo molto importante, perché permette di generare qualsiasi distribuzione/PDF. 

Si parte estraendo un valore $x_i$ in maniera uniforme in $[a,b]$, che è l'intervallo in cui si vuole fare la generazione. I valori sono tutti tenuti, ma sono _pesati_ con $w_i = f(x_i)$. 

La differenza principale riguarda questa volta come è calcolato l'errore sulla distribuzione ottenuta. Questo non può essere trattato come errore poissoniano (ovvero errore di un conteggio semplice). 

:::{admonition} Come trattiamo allora l'errore in questo caso? 
:class: tip
Suppongo che un dato bin dell'istogramma sia composto da N valori, pesati con $w_i$. Allora l'altezza del bin e la sua varianza saranno 

$$
    S = \sum_{i=1}^N w_i \qquad \sigma_S^2 = \sum_{i=1}^N w_i^2.
$$

È utile riscrivere la varianza di S in termini di fluttuazioni del numero di eventi e del peso degli eventi, per osservare alcune caratteristiche di questa descrizione. 

La deviazione standard dei pesi è data come 

$$
    \sigma_w^2 = \langle w^2\rangle - \langle w\rangle^2,
$$ 

dove si ha che $\langle w^2\rangle = \sigma_S^2 /N$ e $\langle w\rangle= S/N$. 

Allora avremo che 

$$
    \frac{\sigma_S^2}{S^2} = \langle w\rangle^2 \frac{N}{S^2} = \frac{N}{S^2} \left( \sigma_w^2 + \langle w\rangle^2 \right) = \frac{N}{S^2} \frac{S^2}{N^2} \left( \frac{\sigma_w^2}{\langle w\rangle^2} + 1 \right),
$$

ovvero che 

$$
    \frac{\sigma_S}{S} = \frac{1}{\sqrt{N}}\sqrt{\frac{\sigma_w^2}{\langle w\rangle^2} + 1}.
$$

Questa riscrittura, che non cambia in alcun modo il risultato che infatti è semplicemente una riscrittura della relazione precedente, rende esplicito il contributo del numero di conteggi e della varianza dei pesi. Se infatti si considerano tutti pesi uguali (e non necessariamenti uguali ad 1) si ottiene che $\sigma_w = 0$, per cui si ottiene che la deviazione standard dei dati è indipendente dal peso che questi hanno. 
:::

