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

##### Teorema della probabilità congiunta

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
    P_C^{(M-3)} (i_{M-3} | i_{M-4}) 
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
    \sum_{i=1}^{\mathcal S} p_i^{(0)} = 1 \quad sum_{j=1}^{\mathcal S} \pi_{ij} = 1, 
$$

corrispondenti alla normalizzazione delle probabilità e anche delle probabilità di transizione. 
:::

























