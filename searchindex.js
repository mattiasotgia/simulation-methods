Search.setIndex({"docnames": ["content/index", "content/lectures/importance_sampling", "content/lectures/introduction", "content/zbiblio"], "filenames": ["content/index.md", "content/lectures/importance_sampling.md", "content/lectures/introduction.ipynb", "content/zbiblio.md"], "titles": ["Metodi di Simulazione applicati alla Fisica", "<span class=\"section-number\">2. </span>Importance sampling (Markov Chain Monte Carlo)", "<span class=\"section-number\">1. </span>Introduzione", "Riferimenti bibliografici"], "terms": {"note": 0, "e": [0, 1, 2], "appunti": [0, 1], "per": [0, 1, 2], "presso": 0, "dipartimento": 0, "universit\u00e0": 0, "degli": [0, 1, 2], "studi": 0, "genova": 0, "tenuto": 0, "dai": [0, 1], "professori": 0, "fabrizio": 0, "parodi": 0, "riccardo": 0, "ferrando": 0, "secondo": [0, 1], "semestr": 0, "dell": [0, 1, 2], "anno": 0, "accademico": 0, "2023": 0, "2024": [0, 1], "laurea": 0, "magistral": 0, "la": [0, 1, 2], "con": [0, 2], "i": [0, 1, 2], "dettagli": 0, "corsi": 0, "unig": 0, "off": [0, 1], "f": [0, 2], "ins": 0, "67330": 0, "si": [0, 1, 2], "prefigg": 0, "fornir": 0, "le": [0, 1, 2], "conoscenz": 0, "base": 0, "sull": 0, "tecnich": 0, "basat": 0, "sul": 0, "metodo": [0, 1, 2], "mont": 0, "carlo": 0, "applicarl": 0, "della": [0, 2], "materia": [0, 1, 2], "interazioni": [0, 2], "fondamentali": 0, "acquisiranno": 0, "competenz": 0, "caten": 0, "markov": 0, "ed": [0, 2], "particolar": [0, 1, 2], "l": [0, 1], "algoritmo": [0, 2], "metropoli": [0, 2], "ii": 0, "transizioni": 0, "fase": 0, "nel": [0, 1, 2], "ga": [0, 1], "reticolar": [0, 1], "iii": 0, "tempo": [0, 2], "continuo": 0, "all": [0, 1, 2], "equilibrio": [0, 1, 2], "fuori": 0, "iv": 0, "crescita": [0, 2], "aggregati": 0, "frattali": 0, "trasporto": 0, "particel": [0, 1, 2], "nella": [0, 1], "interazion": [0, 1, 2], "decadimento": 0, "spazio": [0, 1], "fasi": [0, 1], "parametrica": 0, "un": [0, 1, 2], "rivelator": [0, 2], "esperimento": 0, "composto": [0, 2], "da": [0, 2], "pi\u00f9": [0, 1], "rivelatori": 0, "introduzion": 0, "import": [0, 2], "sampl": [0, 2], "chain": 0, "mattia": 0, "sotgia": 0, "mattiasotgia1": 0, "twitter": 0, "23": 1, "febbraio": 1, "vogliamo": 1, "calcolar": [1, 2], "pratica": 1, "fenomeni": 1, "danneggiamento": 1, "reattori": 1, "nucleari": 1, "nasc": [1, 2], "questo": 1, "ambito": [1, 2], "lo": 1, "alamo": 1, "ovvero": [1, 2], "sullo": 1, "studio": 1, "fission": 1, "nuclear": [1, 2], "il": [1, 2], "problema": 1, "principal": [1, 2], "era": [1, 2], "quello": 1, "risolver": 1, "integrali": 1, "domini": 1, "dimensionalit\u00e0": 1, "molto": [1, 2], "elevata": 1, "In": [1, 2], "questi": [1, 2], "casi": [1, 2], "algoritmi": [1, 2], "come": 1, "simpson": 1, "o": [1, 2], "dei": [1, 2], "trapezi": 1, "fallisc": 1, "rapidament": 1, "infatti": [1, 2], "data": [1, 2], "una": [1, 2], "dimension": 1, "d": [1, 2], "campionando": 1, "10": 1, "punti": 1, "8": 1, "interni": 1, "al": [1, 2], "dominio": 1, "ha": [1, 2], "che": [1, 2], "integrazion": 1, "corrispond": 1, "percentual": 1, "left": [1, 2], "frac": [1, 2], "right": [1, 2], "0": [1, 2], "quad": [1, 2], "text": 1, "grandi": 1, "vediamo": 1, "\u00e8": [1, 2], "necessario": 1, "trovar": 1, "soluzion": [1, 2], "alternativa": 1, "permetta": 1, "non": [1, 2], "aver": [1, 2], "zeri": 1, "divergenz": 1, "chiamiamo": 1, "dello": [1, 2], "bf": 1, "x": 1, "omega": 1, "ci": [1, 2], "troviamo": 1, "uno": [1, 2], "considerando": [1, 2], "funzion": [1, 2], "su": 1, "allora": 1, "potremo": 1, "calcolarn": 1, "valor": [1, 2], "int_": 1, "A": 1, "consideriamo": 1, "p": [1, 2], "tale": [1, 2], "questa": [1, 2], "normalizzata": 1, "1": [1, 2], "definita": 1, "positiva": 1, "ridefinisco": 1, "integral": 1, "precedent": [1, 2], "g": 1, "mi": 1, "accorgo": 1, "facilment": 1, "voler": 1, "medio": 1, "distribuzion": [1, 2], "idea": [1, 2], "vien": 1, "ment": 1, "campionar": 1, "propozionalment": 1, "campiono": 1, "m": 1, "volt": 1, "x_i": [1, 2], "2": [1, 2], "dot": 1, "se": [1, 2], "ho": 1, "campionato": 1, "modo": [1, 2], "corretto": 1, "ottengo": 1, "simeq": 1, "frac1m": 1, "sum_": [1, 2], "obiettivo": 1, "identificar": 1, "adeguata": 1, "ad": [1, 2], "ottimizzar": 1, "render": 1, "vera": 1, "possibil": 1, "uguaglianza": 1, "prima": 1, "anch": [1, 2], "semplic": 1, "scelta": 1, "quella": 1, "considerar": 1, "uniform": [1, 2], "otteniamo": 1, "vm": 1, "pu\u00f2": [1, 2], "esser": [1, 2], "migliorato": 1, "abbiamo": [1, 2], "innanzitutto": 1, "compreso": 1, "sar\u00e0": 1, "molt": [1, 2], "tratteremo": 1, "meccanica": 1, "lagrangiana": 1, "hamiltoniana": 1, "numero": [1, 2], "enorm": 1, "dimensioni": 1, "dovremo": 1, "realt\u00e0": 1, "sono": [1, 2], "casuali": 1, "funzioni": [1, 2], "fornit": 1, "dalla": [1, 2], "fisica": [1, 2], "vedremo": 1, "prend": 1, "nome": [1, 2], "teoria": 1, "suo": 1, "formalismo": 1, "semi": 1, "completo": 1, "rinunceremo": 1, "per\u00f2": [1, 2], "spazi": 1, "infin": 1, "avremmo": 1, "discreto": 1, "finito": 1, "senso": 1, "perch\u00e9": [1, 2], "calcolo": [1, 2], "numerico": 1, "comunqu": [1, 2], "insiem": 1, "mathcal": 1, "": [1, 2], "stati": 1, "indic": 1, "compieremo": 1, "catena": 1, "estrazioni": 1, "ogni": 1, "estrazion": 1, "corrisponder": 1, "sceglier": 1, "del": [1, 2], "sistema": 1, "tratta": 1, "re": 1, "immission": 1, "sommeranno": 1, "total": 1, "totali": 1, "senza": [1, 2], "dimostrarlo": 1, "formalment": 1, "intuitivo": 1, "capir": 1, "quest": 1, "indipendenti": 1, "scorrel": 1, "equiprobabili": 1, "nonostant": 1, "siano": 1, "noi": 1, "serv": 1, "correlata": 1, "sarebb": 1, "troppo": 1, "motivi": 1, "utilizzazion": 1, "correlazion": 1, "sia": [1, 2], "estrarremo": 1, "i_0": 1, "i_1": 1, "i_m": 1, "descrizion": [1, 2], "dice": 1, "esist": [1, 2], "valori": [1, 2], "assenza": 1, "correlazioni": 1, "avremo": [1, 2], "i_2": 1, "invec": 1, "voglio": 1, "ottener": [1, 2], "poter": [1, 2], "cui": [1, 2], "condizionata": 1, "risultati": 1, "precedenti": [1, 2], "vorrei": 1, "risultato": [1, 2], "inlfuenzato": 1, "tutt": [1, 2], "condizioni": 1, "dall": 1, "inizio": 1, "tempi": 1, "solo": 1, "condizion": 1, "immediatament": 1, "incompleti": 1, "esprimer": 1, "termini": [1, 2], "p_c": 1, "i_": 1, "iterando": 1, "3": 1, "inoltr": 1, "stessa": 1, "forma": 1, "indipendentement": 1, "chiameremo": 1, "omogenea": 1, "caso": 1, "ponendoci": 1, "termico": 1, "punto": 1, "considereremo": 1, "estrarr": 1, "stato": [1, 2], "p_i": 1, "j": 1, "pi_": 1, "ij": 1, "nuova": 1, "osserviamo": 1, "vettor": 1, "componenti": 1, "mentr": [1, 2], "matric": 1, "time": 1, "valgono": 1, "corrispondenti": 1, "alla": [1, 2], "normalizzazion": 1, "transizion": 1, "scrittura": 1, "p_": 1, "pi": 1, "_": 1, "semplifica": 1, "notevolment": 1, "chiedo": 1, "qual": [1, 2], "p_j": 1, "sum_i": 1, "matricial": 1, "underlin": 1, "vec": 1, "eq": 1, "lim_": 1, "infti": 1, "quanto": 1, "detto": [1, 2], "precedenza": 1, "saranno": [1, 2], "normalizz": 1, "questometodo": 1, "preved": 1, "numeri": 1, "fare": [1, 2], "raggiung": 1, "vorremo": 1, "boltzmann": 1, "vorremmo": 1, "partir": 1, "partenza": 1, "vuol": [1, 2], "dire": 1, "appropriament": 1, "andar": 1, "macchina": 1, "successivi": 1, "estratti": [1, 2], "avvenga": 1, "dovr\u00f2": 1, "costruir": 1, "adeguatament": 1, "efficenza": [1, 2], "massima": [1, 2], "possibilment": 1, "convergenza": 1, "rapida": 1, "limit": 1, "possiibli": 1, "soluzioni": 1, "possibili": 1, "soddisfano": 1, "effettivament": 1, "infinit": 1, "scegliern": 1, "rispetti": 1, "vincoli": 1, "richiesti": 1, "dia": 1, "tra": 1, "quell": 1, "disposizion": 1, "costruiamo": 1, "converger\u00e0": 1, "richiesta": 1, "soddisfa": 1, "richiest": 1, "ergodicit\u00e0": 1, "capacit\u00e0": 1, "connetter": 1, "diversi": 1, "passi": 1, "raggiunger": 1, "aperiodicit\u00e0": 1, "va": 1, "intrappolar": 1, "cicli": 1, "convergono": 1, "date": 1, "succed": 1, "ergodica": 1, "aperiodica": 1, "ci\u00f2": 1, "cos\u00ec": 1, "impli": 1, "q": 1, "ora": 1, "introdurr": 1, "bilancio": 1, "dettagliato": 1, "implica": 1, "ji": 1, "conservazion": 1, "flusso": 1, "equivalent": 1, "possiamom": 1, "somm": 1, "differenziamo": 1, "elementi": 1, "diagonali": 1, "gli": 1, "diagon": 1, "alpha_": 1, "dove": [1, 2], "simmetrica": 1, "geq": 1, "k": 1, "neq": 1, "ik": 1, "sufficentement": [1, 2], "comporta": 1, "alcun": [1, 2], "sottigliezz": 1, "consideraimo": 1, "n": [1, 2], "spin": 1, "avenso": 1, "frac12": 1, "configurazion": 1, "singola": 1, "s_i": 1, "configurazioni": 1, "dato": [1, 2], "langl": [1, 2], "rangl": [1, 2], "frac1z": 1, "sum": 1, "beta": 1, "h": [1, 2], "addendi": 1, "significativi": 1, "somma": 1, "aspettazion": 1, "energia": 1, "bassa": 1, "alta": 1, "pesati": [1, 2], "energi": 1, "esponent": 1, "sembra": 1, "perfetto": 1, "studiaremo": 1, "consideriammo": 1, "arrai": 1, "siti": 1, "ell": 1, "considero": 1, "occupati": 1, "fissata": 1, "inizial": 1, "passar": [1, 2], "successiva": 1, "regola": 1, "verr\u00e0": 1, "rappresentata": 1, "definiamo": 1, "temepratura": 1, "t": [1, 2], "interagiscono": 1, "intensit\u00e0": 1, "propto": 1, "fatta": 1, "due": 1, "prime": 1, "vicin": 1, "possano": 1, "scambiar": 1, "particella": 1, "coppi": 1, "primi": 1, "vicini": 1, "contar": 1, "devo": 1, "contorno": 1, "spesso": 1, "saaranno": 1, "fissat": 1, "periodich": 1, "utilizzato": 2, "sacco": 2, "branch": 2, "strumento": 2, "negli": 2, "anni": 2, "40": 2, "continua": 2, "evoluzion": 2, "origini": 2, "risalgono": 2, "met\u00e0": 2, "nell": 2, "progetto": 2, "manhattan": 2, "formalizzatori": 2, "enrico": 2, "fermi": 2, "john": 2, "von": 2, "neumann": 2, "stanislaw": 2, "marcin": 2, "ulam": 2, "fu": 2, "inventato": 2, "seguito": 2, "nichola": 2, "constantin": 2, "riferimento": 2, "nota": 2, "tradizion": 2, "nei": 2, "giochi": 2, "azzardo": 2, "omonimo": 2, "sud": 2, "francia": 2, "roulett": 2, "pare": 2, "venn": 2, "convalescent": 2, "malattia": 2, "giocava": 2, "solitario": 2, "cart": 2, "chies": 2, "foss": 2, "probabilit\u00e0": 2, "completarlo": 2, "correttament": 2, "regol": 2, "rendevano": 2, "difficil": 2, "tuttavia": 2, "pens\u00f2": 2, "sarebbero": 2, "potut": 2, "simular": 2, "comput": 2, "disposizioni": 2, "mazzo": 2, "controllar": 2, "quanti": 2, "riusciva": 2, "completar": 2, "poteva": 2, "vincer": 2, "empirico": 2, "ne": 2, "parl\u00f2": 2, "compres": 2, "potenzial": 2, "intuizion": 2, "svilupp": 2, "primo": 2, "permetteva": 2, "generar": 2, "sopratutto": 2, "teorica": 2, "parla": 2, "propriet\u00e0": 2, "sistemi": 2, "bodi": 2, "eccetto": 2, "armonich": 2, "praticament": 2, "impossibil": 2, "analitica": 2, "quindi": 2, "obbligatorio": 2, "utilizzar": 2, "strumenti": 2, "numerici": 2, "processi": 2, "stocastici": 2, "rndm": 2, "walk": 2, "cristal": 2, "atla": 2, "event": 2, "displai": 2, "mathrm": 2, "bar": 2, "sub": 2, "physic": 2, "simul": 2, "initi": 2, "process": 2, "which": 2, "out": 2, "state": 2, "known": 2, "experiment": 2, "simulazion": 2, "produzion": 2, "eventi": 2, "identici": 2, "quelli": 2, "veri": 2, "crear": 2, "maggior": 2, "statistica": 2, "progettar": 2, "ottiimizar": 2, "central": 2, "necessit\u00e0": 2, "casual": 2, "general": 2, "generator": 2, "deve": 2, "convensionali": 2, "altr": 2, "necessari": 2, "rispettano": 2, "nessun": 2, "pattern": 2, "apparent": 2, "x_n": 2, "scorrelato": 2, "x_": 2, "periodo": 2, "ripetizion": 2, "estremament": 2, "lungo": 2, "veloc": 2, "utilizzo": 2, "particl": 2, "lhc": 2, "assumiamo": 2, "adesso": 2, "posseder": 2, "eta": 2, "int_0": 2, "interessa": 2, "distribuzioni": 2, "continu": 2, "gi\u00e0": 2, "visto": 2, "alcuni": 2, "inversion": 2, "integram\u00ecbili": 2, "garantisc": 2, "reiezion": 2, "valido": 2, "qualsiasi": 2, "dipend": 2, "sta": 2, "prevedibil": 2, "priori": 2, "ma": 2, "model": 2, "depend": 2, "rispetto": 2, "miglior": 2, "resta": 2, "omdel": 2, "ancora": 2, "important": 2, "permett": 2, "pdf": 2, "part": 2, "estraendo": 2, "maniera": 2, "b": 2, "intervallo": 2, "tutti": 2, "tenuti": 2, "w_i": 2, "differenza": 2, "riguarda": 2, "volta": 2, "calcolato": 2, "sulla": 2, "ottenuta": 2, "trattato": 2, "poissoniano": 2, "conteggio": 2, "suppongo": 2, "bin": 2, "istogramma": 2, "altezza": 2, "sua": 2, "varianza": 2, "qquad": 2, "sigma_": 2, "util": 2, "riscriver": 2, "fluttuazioni": 2, "peso": 2, "osservar": 2, "caratteristich": 2, "deviazion": 2, "standard": 2, "pesi": 2, "sigma_w": 2, "w": 2, "sqrt": 2, "riscrittura": 2, "cambia": 2, "semplicement": 2, "relazion": 2, "rend": 2, "esplicito": 2, "contributo": 2, "conteggi": 2, "considerano": 2, "uguali": 2, "necessariamenti": 2, "ottien": 2, "dati": 2, "indipendent": 2, "dal": 2, "hanno": 2}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"metodi": [0, 2], "di": [0, 1, 2], "simulazion": 0, "applicati": 0, "alla": 0, "fisica": 0, "il": 0, "corso": 0, "breve": 0, "pagina": 0, "del": 0, "contenuti": 0, "lezioni": 0, "riferimenti": [0, 3], "bibliografici": [0, 3], "contatti": 0, "esterni": 0, "import": 1, "sampl": 1, "markov": 1, "chain": 1, "mont": [1, 2], "carlo": [1, 2], "caten": 1, "teorema": 1, "della": 1, "probabilit\u00e0": 1, "congiunta": 1, "da": 1, "sistemar": 1, "cambio": 1, "notazion": 1, "con": 1, "algoritmo": 1, "metropoli": 1, "modello": 1, "ma": 1, "quindi": 1, "cosa": 1, "andiamo": 1, "studiar": 1, "introduzion": 2, "generazion": 2, "mc": 2, "numeri": 2, "casuali": 2, "campionamento": 2, "semplic": 2, "come": 2, "trattiamo": 2, "allora": 2, "l": 2, "error": 2, "questo": 2, "caso": 2}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 60}, "alltitles": {"Metodi di Simulazione applicati alla Fisica": [[0, "metodi-di-simulazione-applicati-alla-fisica"]], "Il corso in breve": [[0, "il-corso-in-breve"]], "Pagina del corso": [[0, null]], "Contenuti": [[0, "contenuti"]], "Lezioni": [[0, null]], "Riferimenti bibliografici": [[0, null], [3, "riferimenti-bibliografici"]], "Contatti esterni": [[0, null]], "Importance sampling (Markov Chain Monte Carlo)": [[1, "importance-sampling-markov-chain-monte-carlo"]], "Catene di Markov": [[1, "catene-di-markov"]], "": [[1, null], [2, null]], "Teorema della probabilit\u00e0 congiunta": [[1, "teorema-della-probabilita-congiunta"]], "Da sistemare": [[1, null]], "Cambio di notazione": [[1, null]], "Monte Carlo con importance sampling": [[1, "monte-carlo-con-importance-sampling"]], "Algoritmo di Metropolis": [[1, "algoritmo-di-metropolis"]], "Modello": [[1, "modello"]], "Ma quindi cosa andiamo a studiare?": [[1, "ma-quindi-cosa-andiamo-a-studiare"]], "Introduzione": [[2, "introduzione"]], "Metodi di generazione Monte Carlo (MC)": [[2, "metodi-di-generazione-monte-carlo-mc"]], "Numeri casuali": [[2, "numeri-casuali"]], "Campionamento semplice": [[2, "campionamento-semplice"]], "Come trattiamo allora l\u2019errore in questo caso?": [[2, null]]}, "indexentries": {}})