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

Chiameremmo questo caso il caso di __campionamento semplice__, e chiameremo invece ogni altro caso __importance sampling__, dove la distribuzione non è scelta a caso ma sarà scelta in funzione del modello di $A({\bf x})$. 

## Catene di Markov


