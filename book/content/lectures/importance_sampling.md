# Importance sampling (Markov Chain Monte Carlo)

(23 febbraio 2024)

Vogliamo calcolare nella pratica fenomeni di danneggiamento in reattori nucleari. Nasce questo in ambito _Los Alamos_, ovvero sullo studio della fissione nucleare. Il problema principale era quello di risolvere integrali in domini di dimensionalità molto elevata. In questi casi, algoritmi come _simpson_ o il metodo dei trapezi fallisce rapidamente. Infatti data una dimensione $D$, campionando 10 punti (8 punti interni al dominio) si ha che il dominio di integrazione corrisponde in percentuale a 

$$
    \left(\frac{8}{10}\right)^D \to 0 \quad \text{per $D$ grandi.}
$$

Quindi vediamo che è necessario trovare una soluzione alternativa, che permetta di non avere zeri o divergenze. 

Chiamiamo i punti dello spazio $\bf x$, il dominio di integrazione $\Omega$, e ci troviamo in uno spazio di dimensione $D$. 

Considerando una funzione su questo dominio, allora potremo calcolarne il valore come 

$$
    \int_\Omega A({\bf x}) d {\bf x}
$$

Consideriamo una probabilità $p({\bf x})$, tale che questa è normalizzata a 1 e definita positiva. ridefinisco l'integrale precedente come 

$$
    I = \int_\Omega \frac{A({\bf x})}{p({\bf x})} p({\bf x}) d {\bf x} = \int_\Omega G({\bf x}) p({\bf x}) d {\bf x}.
$$

Mi accorgo facilmente che questo corrisponde a voler calcolare il valor medio di $G$ secondo la distribuzione di probabilità $p$. L'idea che allora mi viene in mente è quello di campionare lo spazio delle $\bf x$ propozionalmente a $p$, lo campiono $M$ volte, ${\bf x_i}$, $i=1, 2,\dots,M$. Se ho campionato in modo corretto ottengo che 

$$
    I\simeq \frac1M \sum_{i=1}^M G({\bf x_i}) = \frac1M \sum_{i=1}^M \frac{A({\bf x_i})}{p({\bf x_i})}
$$

Allora l'obiettivo principale è quello di identificare una $p$ adeguata ad ottimizzare e rendere più vera possibile questa uguaglianza. La prima (e anche più semplice) scelta che mi viene è quella di considerare una distribuzione uniforme. Otteniamo allora che 

$$
    I\simeq \frac VM \sum_{i=1}^M A({\bf x_i}).
$$

:::{admonition} 
:class: attention
Chiameremmo questo caso il caso di __campionamento semplice__, e chiameremo invece ogni altro caso __importance sampling__, dove la distribuzione non è scelta a caso ma sarà scelta in funzione del modello di $A({\bf x})$. 
:::

## Catene di Markov

Il modello precedente può essere molto migliorato. Abbiamo innanzitutto compreso che 
 1. Lo spazio sarà a molte dimensione (tratteremo lo spazio delle fasi secondo la meccanica lagrangiana e hamiltoniana, quindi con un numero enorme di dimensioni)
 2. Dovremo avere probabilità che in realtà non sono casuali, ma funzioni $p({\bf x})$ fornite dalla fisica.

Vedremo quindi quella che prende il nome dei teoria delle catene di Markov. La tratteremo nel suo formalismo (semi) completo. Rinunceremo però a considerare spazi infiniti, ma avremmo uno spazio delle fasi __discreto__ e __finito__. Questo ha molto senso perché il calcolo numerico è comunque un calcolo discreto e finito. 

Consideriamo un insieme di $\mathcal S$ stati, e l'indice sarà $i = 1, \dots, \mathcal S$. Su questo spazio di stati compieremo una __catena di estrazioni__. Ogni estrazione corrispondere a scegliere uno degli stati del sistema, e si tratta di estrazioni con re-immissione. Consideriamo $0, \dots, M$ estrazioni, che sommeranno ad un totale di $M+1$ estrazioni totali. 

:::{margin}
Senza dimostrarlo formalmente, è intuitivo capire che queste estrazioni sono in realtà indipendenti (ovvero scorrelate) e equiprobabili.
:::

Nonostante le Markov Chains siano scorrelate, quello che in realtà a noi serve è una catena correlata (sarebbe troppo semplice), e vogliamo che per motivi di utilizzazione pratica questa correlazione sia la più semplice possibile. 

Nella catena di estrazioni estrarremo 

$$
    i_0, i_1, \dots, i_M.
$$

La descrizione di questa estrazione mi dice che esiste la probabilità __congiunta__ a $M+1$ valori 

$$
    p^{(M)} (i_0, i_1, \dots, i_M).
$$

In assenza di correlazioni avremo che 

$$
    p^{(M)} (i_0, i_1, \dots, i_M) = p^{(1)} (i_0) p^{(1)} (i_1) p^{(1)} (i_2) \dots p^{(1)} (i_M)
$$

Quello che invece voglio ottenere è di poter calcolare delle catene in cui l'estrazione è condizionata dai risultati precedenti. Non vorrei però che il risultato sia inlfuenzato da tutte le condizioni precedenti dall'inizio dei tempi, ma solo dalla condizione immediatamente precedente. 

### Teorema della probabilità congiunta

:::{admonition} Da sistemare
:class: warning
Appunti incompleti
:::

È possibile esprimere la probabilità congiunta in termini di una probabilità condizionata 

$$
    p^{(M)} (i_0, i_1, \dots, i_M) = P_C^{(M)} (i_M | i_0, i_1, \dots, i_{M-1}) p^{(M-1)} (i_0, i_1, \dots, i_{M-1})
$$

Iterando allora ottengo che 

$$
    p^{(M)} (i_0, i_1, \dots, i_M) = 
    P_C^{(M)} (i_M | i_{M-1}) 
    P_C^{(M-1)} (i_{M-1} | i_{M-2}) 
    P_C^{(M-2)} (i_{M-2} | i_{M-3}) 
    \dots 
    P_C^{(1)} (i_{1} | i_{0}) 
    p^{(0)} (i_0)
$$

Se inoltre $P_C^{(m)}$ ha la stessa forma indipendentemente da $m$, allora chiameremo questa catena una catena di Markov omogenea. Questo sarà il caso principale che tratteremo, ponendoci all'equilibrio termico. 


:::{admonition} Cambio di notazione
:class: warning
Da questo punto considereremo un cambio di notazione. 

Chiamiamo la probabilità di estrarre lo stato $i$ all'inizio come 

$$
    p^{(0)}(i_0) \to p_i ^{(0)}.
$$

Inoltre avremo che 

$$
    P_C^{(M)} (j | i) \to \pi_{ij}
$$

Con questa nuova notazione osserviamo che $p_i^{(0)}$ è un vettore con $\mathcal S$ componenti, mentre $\pi_{ij}$ è una matrice di dimensione $\mathcal S \times \mathcal S$. Valgono in questa nuova notazione le condizioni 


$$
    \sum_{i=1}^{\mathcal S} p_i^{(0)} = 1 \quad \sum_{j=1}^{\mathcal S} \pi_{ij} = 1, 
$$

corrispondenti alla normalizzazione delle probabilità e anche delle probabilità di transizione. 
:::

Avremo che con questa nuova notazione la scrittura 

$$
    p^{(1)} (i_0, i_1) = P_C^{(1)} (i_1 | i_0) p^{(0)} (i_0) \to p_{ij}^{(1)} = p_i^{(0)} \pi _{ij},
$$

ovvero si semplifica notevolmente. MA allora se mi chiedo quale è la probabilità di estrarre lo stato $j$ all'estrazione 1, questa sarà data da 

$$ 
    p_j^{(1)} = \sum_i p_i^{(0)} \pi _{ij},
$$ 

ovvero in notazione matriciale ${\bf p}^{(1)} = {\bf p}^{(0)} \underline \pi$. 

Iterando allora ottengo che $\vec p^{(m)} = \vec p^{(0)} \underline \pi^{m}$, che mi dice allora che 

$$ 
    \vec p^\text{eq} = \lim_{m\to\infty} \vec p^{(0)} \underline \pi^{m}.
$$


Osserviamo che per quanto detto in precedenza, allora tutte le probabilità saranno normalizzate ad uno. 

### Monte Carlo con _importance sampling_

Questometodo prevede di estrarre numeri per fare si che la distribuzione che si raggiunge sia la distribuzione di equilibrio, che vorremo essere in particolare per la fisica della materia distribuzione di Boltzmann. Questa la chiameremo $\vec p^\text{eq}$. Vorremmo che 

$$ 
    \vec p^\text{eq} = \lim_{m\to\infty} \vec p^{(0)} \underline \pi^{m}.
$$

a partire da ogni probabilità $p^{(0)}$ di partenza. Questo vuol dire 
 1. Scegliere appropriamente una $p^{(0)}$
 2. Fare andare la macchina, in modo che campionando stati successivi siano estratti secondo la probabilità $p^\text{eq}$.

Perché questo avvenga dovrò costruire la matrice $\pi$ adeguatamente per avere l'efficenza massima (ovvero possibilmente per una convergenza rapida del limite). 

Le possiibli soluzioni, ovvero le possibili $\pi$ che soddisfano questo limite, sono effettivamente infinite, e quello che noi vorremo fare è di sceglierne una in particolare che rispetti i vincoli richiesti, ma che dia una soluzione semplice tra quelle a disposizione. 

La $\pi$ che costruiamo convergerà alla richiesta $\vec p^\text{eq}$ se soddisfa richieste di

 1. Ergodicità. Capacità di connettere stati diversi in un numero finito di passi successivi. Nel caso delle $\pi$ vuol dire che a partire da ogni stato di partenza $p^{(0)}$ sarà possibile raggiungere $\vec p^\text{eq}$. 
 2. Aperiodicità. Non si va ad intrappolare in cicli che non convergono.

Date queste condizioni, quello che succede è che 

__Teorema__. Se $\pi$ è _ergodica_ e _aperiodica_, allora esiste $\vec p^*$ tale che 

$$
    \lim_{m\to\infty} \vec p^{(0)} \pi^m = \vec p^*,
$$

per ogni $\vec p^{(0)}$. 

Vogliamo che $\vec p^\text{eq} = \vec p^*$. Per avere ciò è necessario scegliere $\pi$ per cui $\vec p^\text{eq} \pi = \vec p^\text{eq}$. Infatti così ho che 

$$
    \vec p^\text{eq} \pi^m = \vec p^\text{eq} \implies \lim_{m\to\infty}\vec p^\text{eq} \pi = \vec p^\text{eq} \implies \text{Q.E.D.}
$$

Voglio costruire ora la $\pi$. Per fare ciò voglio introdurre il __bilancio dettagliato__. Questo implica una scelta di $\pi$ data da 

$$
    \vec p^\text{eq}_{i} \pi_{ij} = \vec p^\text{eq}_{j} \pi_{ji}.
$$

Questo corrisponde ad una conservazione del flusso. La transizione $i\to j$ è equivalente a $j\to i$.

Per trovare allora la forma di $\pi$ possiamom considerare le somme su $i$, ovvero 

$$
    \sum_i \vec p^\text{eq}_{i} \pi_{ij} = \sum_i \vec p^\text{eq}_{j} \pi_{ji} = \vec p^\text{eq}_{j} \sum_i \pi_{ji} = \vec p^\text{eq}_{j}.
$$

## Algoritmo di Metropolis

Differenziamo tra elementi non diagonali e gli elementi diagonali.

 - Elementi off-diagonal.

   $$
       \vec p^\text{eq}_{j} < \vec p^\text{eq}_{j} \implies \pi_{ij} = \alpha_{ij} \frac{\vec p^\text{eq}_{j}}{\vec p^\text{eq}_{i}}
   $$
   
   Dove $\alpha_{ij}$ è una matrice simmetrica, e 

   $$
       \vec p^\text{eq}_{j} \geq \vec p^\text{eq}_{j} \implies \pi_{ij} = \alpha_{ij}
   $$

 - Elementi diagonal

   $$
       \pi_{ij} = 1 - \sum_{k\neq i} \pi_{ik}
   $$

Questa scelta soddisfa il bilancio dettagliato. 













































