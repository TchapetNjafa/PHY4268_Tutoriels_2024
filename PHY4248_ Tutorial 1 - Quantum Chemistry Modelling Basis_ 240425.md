---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.6
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# PHY4268 Tutorial 1 - Quantum Chemistry Modelling Basis

1. **S. G. Nana Engo**, serge.nana-engo@facsciences-uy1.cm
    * Department of Physics, Faculty of Science, University of Yaounde I
2. **J-P. Tchapet Njafa**, jean-pierre.tchapet-njafa@univ-maroua.cm
    * Department of Physics, Faculty of Science, University of Maroua
       
April 2024


Dans cette série de tutoriels, nous allons aborder les différentes composantes d'un algorithme de la VQE permettant de faire des simulations moléculaires, comme illustrée par la Figure ci-jointe.

![VQE_Flowchart_Tuto1.jpeg](./Graphics/VQE_Flowchart_Tuto1.jpeg)

Comme nous l'avons déjà précédemment vu, le calcul quantique impacte de nombreux secteurs technologiques et pharmacetiques. La Figure ci-dessous pour présente les applications dans le domaine de l'énergie.

![CoverVQE.jpg](./Graphics/CoverVQE.jpg)

A l'issue de cette série de tutoriels, en utilisant les frameworks de calculs de chimie quantique **PySCF** et de cacul quantique **qiskit**, l'étudiant doit être capable 
* d'expliquer certains concepts des méthodes modernes de chimie quantique : Hartree-Fock et la théorie fonctionnelle de la densité (DFT) ;
* d'interpréter les résultats calculés en termes d'orbitales moléculaires et de modèles de structure électronique ;
* d'implémenter l'algorithme de la VQE avec qiskit-nature;
* d'appliquer cet algorithme pour prédire les structures moléculaires.



L'objectif de ce **tutoriel 1** est de fournir les connaissances de base en chimie quantique nécessaire résoudre un problème de chimie avec HF et la DFT (voir les parties marquées en rouge sur la Figure 1 présentant la VQE. De nombreux paramètres peuvent être modifiés, et certains ont un impact plus important que d’autres sur la précision et les ressources computationnelles. Cet ensemble de paramètres comprend le choix de l'état initial, le choix de l'ensemble de bases. De façon spécifique, à la fin de ce tutoriel, l'apprenant doit
* avoir compris les notions théoriques élémentaires de la structure électronique des molécules en seconde quantification ;
* être capable d'utiliser le framework PySCF pour 
    * effectuer le calcul du champ moyen (SCF)
    * extraire
        * l'énergie totale
        * les index et les énergies du HOMO et du LUMO, 
        * le résume des différentes énergies calculées
    * explorer le choix du set de bases.


<!-- #region -->
## Structure électronique des molécules en première et seconde quantifications

### Hamiltonien d'une molécule en première quantification

L'Hamiltonien d'une molécule constituée de noyaux $K$ et d'électrons $N$ est
\begin{equation*}
\mathtt{H}=\mathtt{T}_e+\mathtt{T}_n+\mathtt{U}_{en}+\mathtt{U}_{en}+\mathtt{U}_{nn},
\end{equation*}
où
- $\mathtt{T}_e= -\sum_i\frac{\hbar^2}{2m_e}\nabla^2_i$ est l'énergie cinétique électronique;
- $\mathtt{T}_n=-\sum_I\frac{\hbar^2}{2M_I}\nabla^2_I$ est l'énergie cinétique nucléaire;
- $\mathtt{U}_{en}= -\sum_{i,I}\frac{e^2}{4\pi\epsilon_0}\frac{Z_I}{|\mathbf{r}_i-\mathbf{R}_I|}$ est la répulsion Coulombienne entre les électrons et les noyaux;
- $\mathtt{U}_{ee}=+\frac{1}{2}\sum_{i\neq j}\frac{e^2}{4\pi\epsilon_0}\frac{1}{|\mathbf{r}_i-\mathbf{r}_j|}$ est la répulsion Coulombienne entre les électrons eux-mêmes;
- $\mathtt{U}_{nn}=+\frac{1}{2}\sum_{I\neq J}\frac{e^2}{4\pi\epsilon_0}\frac{Z_IZ_J}{|\mathbf{R}_I-\mathbf{R}_J|}$ est la répulsion Coulombienne entre les noyaux eux-mêmes.

$M_I$, $\mathbf{R}_I$ et $Z_I$ ($I=1,2$) désignent la masse, la position et le numéro atomique du *I*-ème noyau, et $\mathbf{r}_i$ est la position du *i*-ème électron. Par souci de concision, nous travaillons en unités atomiques, où l'unité de longueur est $a_0 = 1$ Bohr ($0.529167~\times~10^{-10}$ m), l'unité de masse est la masse électronique $m_e$ et l'unité d'énergie est 1 Hartree ($1~\textrm{Hartree}~= ~{e^2}/{4\pi\epsilon_0 a_0}~=~27.2113~\textrm{eV}$). Dénotant $M_I'=M_I/m_e$, les différents termes de l'Hamiltonien moléculaire en unités atomiques s'écrit
\begin{equation*}
\begin{split}
\mathtt{T}_e&=-\sum_i\frac{\nabla^2_i}{2}, \qquad\mathtt{T}_n= -\sum_I\frac{\nabla^2_I}{2M_I'}, 
 \qquad\mathtt{U}_{en}= - \sum_{i,I}\frac{Z_I}{|\mathbf{r}_i-\mathbf{R}_I|},\\
\mathtt{U}_{ee}&=+\frac12\sum_{i\neq j}\frac{1}{|\mathbf{r}_i-\mathbf{r}_j|}
=+\sum_{i,j>i}\frac{1}{|\mathbf{r}_i-\mathbf{r}_j|},\\
\mathtt{U}_{nn}&=+\frac{1}{2}\sum_{I\neq J}\frac{Z_IZ_J}{|\mathbf{R}_I-\mathbf{R}_J|}
=+\sum_{I,J>I}\frac{Z_IZ_J}{|\mathbf{R}_I-\mathbf{R}_J|}.
\end{split}
\end{equation*}

### Approximation de Born-Oppenheimer

Nous nous intéressons principalement à la structure *électronique* de la molécule. Comme $M'\sim10^3$, nous appliquons l'**approximation de Born-Oppenheimer**, en traitant les noyaux comme des charges ponctuelles classiques. En conséquence, le terme d'énergie cinétique nucléaire $T_n$ est négligé, $\mathtt{U}_{nn}(\mathbf{R})$ est constant et l'Hamiltonien moléculaire est approché par un Hamiltonien électronique paramétré par $\mathbf{R}$,
\begin{equation*}
\begin{split}
\mathtt{H}(\mathbf{R}) &= \mathtt{U}_{nn}+\underset{\text{One-Electron Operators}}{\mathtt{T}_e + \mathtt{U}_{en}(\mathbf{R})} + \underset{\text{Two-Electron Operator}}{\mathtt{U}{ee}}\\
&= \mathtt{U}_{nn}+\underset{\text{One-Electron Operator}}{\sum_i h(i)} + \underset{\text{Two-Electron Operator}}{\sum_{i\ne j}v(i, j)}.
\end{split}
\end{equation*}

### Énergie potentiel de surface

En vertu de l'équation de $\mathtt{H}(\mathbf{R})$ préceédente, l'équation de Schrödinger indépendante du temps non relativiste écrit
\begin{equation*}
\mathtt{H}(\mathbf{R})\psi_i(\mathbf{R},\mathbf{r})=E_i(\mathbf{R})\psi_i(\mathbf{R},\mathbf{r}),
\end{equation*}
où $\psi_i(\mathbf{R},\mathbf{r})$ désigne l'un des états propres et $E_i(\mathbf{R})$ correspond aux surfaces d'énergie potentielle moléculaire (PES), ce qui est important pour comprendre les processus chimiques tels que la formation et la rupture de liaison. En règle générale, la précision du calcul du PES (Potential Energy Surface) doit être limitée à la *précision chimique* $\epsilon=1 \textrm{kcal}\textrm{mol}^{-1}(\sim 1.59\times 10^{-3}\textrm{Hartrees},\, 43.3\textrm{meV})$.

L'énergie de l'état fondamental est donnée par
\begin{equation}\tag{1} 
E_0 = \frac{\langle\psi_0|\mathtt{H}|\psi_0\rangle}{\langle\psi_0|\psi_0\rangle},
\end{equation}
où $|\psi_0\rangle$ est l'état fondamental du système. Cependant, la dimensionnalité de ce problème croît exponentiellement avec le nombre de degrés de liberté. Pour résoudre ce problème, nous aimerions préparer $|\psi_0\rangle$ sur un calculateur quantique et mesurer directement la valeur moyenne Hamiltonienne (ou $E_0$). Alors, comment pouvons-nous faire cela concrètement?

Un bon point de départ pour résoudre ce problème est la méthode Hartree-Fock (HF), autrement nous allons prendre comme fonction d'état de réféfrence la fonction d'état HF. Cette méthode rapproche un problème à N corps en N problèmes à un corps où chaque électron évolue dans le champ moyen des autres.  Il est à noter que le dénominateur de l'Equation (1) n'est nécessaire que quand la fonction d'etat HF n'est pas normalisée au départ. Une facon d'améliorer cette fonction d'état HF est d'utiliser une combinaison linéaire des fonction d'état. Ce qui sera fait à travers les bases chimiques.


### Déterminant Slater

Un déterminant de Slater est un produit anti-symétrique d'un ou de plusieurs spin orbitals. Pour un système à $N$-electron, le déterminant de Slater s'écrit
\begin{equation*}
\begin{split}
\psi(\mathbf{r}_0, \mathbf{r}_1, \ldots, \mathbf{r}_{N-1})\equiv |\phi_{M-1},\cdots,  \phi_1, \phi_0 \rangle 
&=   \frac{1}{\sqrt{N!}}
\begin{vmatrix} \phi_0(\mathbf{r}_0) & \phi_1(\mathbf{r}_0) & \cdots & \phi_{M-1}(\mathbf{r}_0) \\
      \phi_0(\mathbf{r}_1) & \phi_1(\mathbf{r}_1) & \cdots & \phi_{M-1}(\mathbf{r}_1) \\
      \vdots & \vdots & \ddots & \vdots \\
      \phi_0(\mathbf{r}_{N-1}) & \phi_1(\mathbf{r}_{N-1}) & \cdots & \phi_{M-1}(\mathbf{r}_{N-1})
\end{vmatrix} \\
&\equiv |M-1,\dots, 1, 0 \rangle,
\end{split}
\end{equation*}
 L'échange des positions de deux électrons est équivalent à l'échange de deux rangées du déterminant de Slater, ce qui modifie le signe de la fonction d'état. Cela fournit la symétrie d'échange correcte pour la fonction d'état fermionique. Bien que le nombre de spin-orbitales considérées, $M$, soit généralement plus élevé que le nombre d'électrons dans le système, $N$, les électrons ne peuvent occuper que $N$ des spin-orbitales dans un déterminant Slater donné. En conséquence, le déterminant de Slater ne contient que les spin-orbitales occupées $N$.
<!-- #endregion -->

### Exemple d'illustration

Pour rendre les correspondances plus claires, nous examinons comment une fonction d'état sera représentée pour un système quantique fictif. On considère les orbitales de spin $|A_\uparrow\rangle,\,|A_\downarrow\rangle,\,|B_\uparrow\rangle,\,|B_\downarrow\rangle$. Nous sommes libres de définir arbitrairement l'état Hartree-Fock de notre système fictif et de choisir les deux électrons dans l'orbitale $|A\rangle$. Nous sommes intéressés par la fonction d'état lorsque la composante $z$ du spin est nulle.

Labellisons chaque orbitale ainsi qu'il suit :

\begin{align*}
 &|A_\uparrow\rangle=|00\rangle=|\mathbf{0}\rangle,  &&|A_\downarrow\rangle=|01\rangle=|\mathbf{1}\rangle,
&&|B_\uparrow\rangle=|10\rangle=|\mathbf{2}\rangle,  &|B_\downarrow\rangle=|11\rangle=|\mathbf{3}\rangle.
\end{align*}

L'état de Hartree-Fock possède les deux électrons dans les orbitales $|A\rangle$. Un état HF incorrectement symétrisé serait donc $|A_\uparrow\rangle_1|A_\downarrow\rangle_2 = |\mathbf{0}\rangle_1|\mathbf{1}\rangle_2$, où les indices indiquent l'électron décrit par chaque orbitale. La fonction d'état HF correctement antisymétrisée serait 

\begin{equation*}
\begin{split}
    |\Psi_{\mathrm{HF}}\rangle =  \frac{1}{\sqrt{2}}
		\begin{vmatrix}
		A_\uparrow(\mathbf{x}_0) & A_\downarrow(\mathbf{x}_0)\\
		A_\uparrow(\mathbf{x}_1) & A_\downarrow(\mathbf{x}_1)\\
		\end{vmatrix}
		&= \frac{1}{\sqrt{2}}  \Big( A_\uparrow(\mathbf{x_0}) A_\downarrow(\mathbf{x_1}) - A_\downarrow(\mathbf{x_0}) A_\uparrow(\mathbf{x_1})  \Big)\\
&=\frac{1}{\sqrt{2}} (|\mathbf{0}\rangle_1|\mathbf{1}\rangle_2 - |\mathbf{1}\rangle_1|\mathbf{0}\rangle_2 ).
\end{split}
\end{equation*}

Si nous considérons maintenant les excitations au-dessus de l'état HF, alors une fonction d'état générale avec $s_z=0$ qui a été correctement antisymétrisée est
donnée par

\begin{equation*}
\begin{aligned}
|\Psi\rangle =& \frac{\alpha}{\sqrt{2}}(|\mathbf{0}\rangle_1|\mathbf{1}\rangle_2 - |\mathbf{1}\rangle_1|\mathbf{0}\rangle_2 ) 
+\frac{\beta}{\sqrt{2}} (|\mathbf{2}\rangle_1|\mathbf{3}\rangle_2 - |\mathbf{3}\rangle_1|\mathbf{2}\rangle_2 ) \\
+& \frac{\gamma}{\sqrt{2}} (|\mathbf{0}\rangle_1|\mathbf{3}\rangle_2 - |\mathbf{3}\rangle_1|\mathbf{0}\rangle_2 ) + \frac{\delta}{\sqrt{2}} (|\mathbf{1}\rangle_1|\mathbf{2}\rangle_2 - |\mathbf{2}\rangle_1|\mathbf{1}\rangle_2 ).
\end{aligned}
\end{equation*}

Comme nous avons $N=2$ électrons et $M=4$ spin-orbitales, nous pouvons voir que nous avons seulement besoin de $N \lceil\mathrm{log_2}(M)\rceil = 2 \times \lceil\mathrm{log_2}(4)\rceil = 4$ qubits pour stocker la fonction d'état.



## Seconde quantification - Construction d'un opérateur Hamiltonien fermionique

### Génération des intégrales moléculaires

L'Hamiltonien  est exprimé dans la base des solutions de la méthode HF, également appelées Orbitales Moléculaires (OM) :
$$
\mathtt{H}_{elec}=\sum_{pq} h_{pq} a^{\dagger}_p a_q + 
\frac{1}{2} \sum_{pqrs} h_{pqrs}  a^{\dagger}_p a^{\dagger}_q a_r  a_s
$$
avec 
* les **intégrales à 1 électron**
$$
h_{pq} = \int \phi^*_p(r) \left( -\frac{1}{2} \nabla^2 - \sum_{I} \frac{Z_I}{|\mathbf{r}_i-\mathbf{R}_I|} \right)   \phi_q(r)d\mathbf{r}
$$
qui décrivent l’énergie cinétique des électrons individuels et leurs interactions Coulombiennes avec les champs électriques du noyau ;
* et **intégrales à 2 électrons**
$$
h_{pqrs} = \int \frac{\phi^*_p(r_1)  \phi^*_q(r_2) \phi_r(r_2)  \phi_s(r_1)}{|\mathbf{r}_1-\mathbf{r}_2|}d\mathbf{r}_1d\mathbf{r}_2,
$$
décrivent les interactions Coulombiennes entre les électrons;
* $a^\dagger_p$ l'**opérateur de création** d'un électron sur le spin-orbite $p$;
* $a_q$ l'**opérateur d'annihilation** d'un électron sur sur le spin-orbite $q$;
* $a^{\dagger}_p a_q $ est l'**opérateur d'excitation**, qui excite un électron du spin-orbite occupé $\phi_q$ vers le spin-orbite inoccupé $\phi_p$.

Les MO ($\phi_u$) peuvent être occupés ou virtuels (inoccupés). Un MO peut contenir 2 électrons. Ces MO sont les solutions de la méthode HF, où chaque l´electron est soumis au champ moyen créé par les autres électrons. Cependant, dans ce qui suit, nous travaillons en fait avec des orbitales de spin qui sont associées à un spin up ($\alpha$) d'électron spin down ($\beta$). Ces deux spins sont également communément désignés par $\alpha$ et $\beta$, respectivement. Ainsi, les orbitales de spin peuvent contenir un électron ou être inoccupées.

Dans le tableau ci-dessous resume les correspondances première-seconde quantification des divers éléments du Hamiltonien moléculaire. 
$$
\begin{array}{|l|l|l|}\hline
\textbf{Approximation BO}  &  \textbf{Combinaison en 2e quantification} &  \textbf{Description}
\\\hline
-\frac12\sum_i\nabla^2_i & -\frac12\sum_{p,q}\langle p|\nabla^2_p|q\rangle a^\dagger_p a_q & \mathtt{T}_e \\\hline
 -\sum_{i,I}\frac{Z_I}{|\mathbf{r}_i-\mathbf{R}_I|} &
 -\sum_{p,q} \langle p|\frac{Z_I}{|\mathbf{r}_p-\mathbf{R}_I|}|q\rangle a^\dagger_p a_q
 & \mathtt{U}_{en}\\\hline
\sum_{i,j>i}\frac{1}{|\mathbf{r}_i-\mathbf{r}_j|} & 
\sum_{p,q,r,s}\langle pq \big|\frac{1}{|\mathbf{r}_i-\mathbf{r}_j|}\big|rs\rangle a^\dagger_p a^\dagger_q a_r a_s & 
 \mathtt{U}_{ee}\\\hline
\end{array}
$$

*Un des avantages apporté par le formalisme de seconde quantification est que la propriété d'anti-symétrie des fonctions d'état pour l'échange de fermions identiques, prise en compte manuellement dans une approche de première quantification en considérant les déterminants de Slater, est automatiquement imposée par le relations d'anti-commutation des opérateurs de création et d'annihilation. Cependant, le prix à payer est de travailler avec un nombre non fixe de particules, qui reste de toute façon une quantité conservée.*


### Représentation nombre d'occupations

En seconde quantification, le déterminant de Slater est représenté par un **vecteur nombre d'occupations**
\begin{equation*}
|f\rangle =|f_{M-1},\dots,f_{i-1},f_i,f_{i+1},\dots, f_0\rangle, 
\end{equation*}
avec
- $f_i=1$ lorsque le spin-orbital $\phi_i$ est occupé et donc présent dans le déterminant de Slater;
- $f_i = 0$ lorsque le $\phi_i$ spin-orbital est vide et donc absent du déterminant de Slater.

- Le vecteur $|f_i\rangle$ est appelé **vecteur nombre d'occupations de l'orbite fermionique i**, parce que
\begin{align*}
& |0\rangle=\begin{pmatrix}1\\0\end{pmatrix}:\text{état vide}, & |1\rangle=\begin{pmatrix}0\\1\end{pmatrix}:\text{état occupé}.
\end{align*}

Dans la deuxième représentation de quantification, au lieu de poser la question *Quel électron est dans quel état?*, nous posons la question *combien de particules y a-t-il dans chaque état?*.

### Opérateurs fermioniques

Les opérateurs fermioniques $a^\dagger_i$ et $a_i$ obéissent aux relations d'anticommutation fermioniques suivantes:
\begin{align*}
&\{a_p,a^\dagger_q\} = a_pa^\dagger_q + a^\dagger_q a_p = \delta_{pq}, 
&\{a_p, a_q \} = \{a^\dagger_p, a^\dagger_q\} = 0 .
\end{align*}
La nature fermionique des électrons implique que les fonctions d'état à plusieurs électrons doivent être antisymétriques par rapport à l'échange de particules. Cela se reflète dans la manière dont les opérateurs de création fermionique et d'annhilation agissent sur les déterminants $|f\rangle$:
\begin{equation*}
\begin{aligned}
&a_i|f_{M-1},\dots,f_{i-1},f_i,f_{i+1},\dots, f_0\rangle = \delta_{f_i,1}
(-1)^{\sum_{m<i}f_m}|f_{M-1},\dots,f_{i-1},f_i\oplus 1,f_{i+1},\dots, f_0\rangle,\\
&a^\dagger_i|f_{M-1},\dots,f_{i-1},f_i,f_{i+1},\dots, f_0\rangle = \delta_{f_i,0}
(-1)^{\sum_{m<i}f_m}|f_{M-1},\dots,f_{i-1},f_i\oplus 1,f_{i+1},\dots, f_0\rangle`,  \\
&a|0\rangle=a^\dagger|1\rangle=0.
\end{aligned}
\end{equation*}
- Le terme de phase $p_i=(-1)^{\sum_{m<i}f_m}$ désigne la parité qu'impose l'anti-symétrie d'échange des fermions.
    - $p_i=-1$ si le nombre d'électrons est impair dans ces orbitales de spin,
    - $p_i= 1$ si le nombre d'électrons est pair dans ces orbitales de spin.

- Le symbole $\oplus$ représente l'addition modulo 2, i.e.,  $1 \oplus 1=0 $ et $0 \oplus 1=1$.

- L'opérateur de création $a^\dagger_i$ met un électron dans une orbitale inoccupée *i* ou dans $|f_i\rangle$.

- L'opérateur d'annihilation $a_i$ supprime un électron dans une orbitale occupée *i* ou en $|f_i\rangle$.

- $a_i^2 =(a^\dagger_i)^2= 0$ pour tout $i$. On ne peut pas créer ou anéantir un fermion dans le même mode deux fois.

Par example,
\begin{align*} 
&a^\dagger_1 |0\rangle_1 = |1\rangle_1,&&a^\dagger_1 |1\rangle_1 =(a^\dagger_i)^2|0\rangle_1= 0,&& a_1 |1\rangle_1 =|0\rangle_1,&a_1 |0\rangle_1 =a_1^2|1\rangle_1=0. 
\end{align*}
Notez qu'ici $a^\dagger_1|1\rangle_1=0$ et $a_1 |0\rangle_1=0$ signifie le vecteur zéro et non $|0\rangle_1$. 

En utilisant également de tels opérateurs, nous pouvons exprimer
\begin{equation*}
|0\rangle |1\rangle |1\rangle |0\rangle = a^\dagger_1a^\dagger_2 |0\rangle^{\otimes 4}.
\end{equation*}

 L'opérateur d'occupation orbitale est donné par
\begin{align*}
&n_i = a^\dagger_i a_i, &n_i |f_{M-1},\dots,f_i,\dots, f_0\rangle= f_i |f_{M-1},\dots,f_i,\dots, f_0\rangle,
\end{align*}
et il compte le nombre d'électrons dans une orbitale donnée, c'est-à-dire
\begin{align*} 
&n_i |0\rangle_i = 0, &n_i |1\rangle_i = |1\rangle_i. 
\end{align*}

### Exercice

In **Exemple d'illustration** above, we have considered a fictitious system described by spin orbitals $|A_\uparrow\rangle,\,|A_\downarrow\rangle,\,|B_\uparrow\rangle,\,|B_\downarrow\rangle$ labelled as
\begin{align*}
 &|A_\uparrow\rangle=|00\rangle=|\mathbf{0}\rangle,  &&|A_\downarrow\rangle=|01\rangle=|\mathbf{1}\rangle,
&&|B_\uparrow\rangle=|10\rangle=|\mathbf{2}\rangle,  &|B_\downarrow\rangle=|11\rangle=|\mathbf{3}\rangle.
\end{align*}
We assume that the Hartree-Fock state for this system **has both electrons occupying the $|A\rangle$ orbitals**. We store the occupations of the spin orbitals $|A_\uparrow\rangle,\,|A_\downarrow\rangle,\,|B_\uparrow\rangle,\,|B_\downarrow\rangle$ which we order as 
$|f_{B_\downarrow},f_{B_\uparrow},f_{A_\downarrow},f_{A_\uparrow}\rangle$, with $f_i=0,1$ (in the fermionic Fock space).

1. Write down, in the fermionic Fock space, the Hartree-Fock state $|\Psi_{\rm HF}\rangle$, that is, the corresponding of the antisymmetrized Slater determinant obtained in **Exemple d'illustration**.
$$
|\Psi_{\rm HF}\rangle = |0011\rangle .
$$
1. Write down, in the fermionic Fock space, the $s_z = 0$ state function.
$$
|\Psi_{\rm HF}\rangle = \alpha|0011\rangle + \beta|1100\rangle + |1001\rangle + |0101\rangle .
$$



## Calcul du champ moyen ou mean-field

Le calcul du champ moyen est généralement la première étape de toutes les approches post Hartree-Fock. Les approches en champ moyen résolvent l'équation de Schrödinger en faisant la moyenne du terme de répulsion électron-électron de manière auto-cohérente ou SCF. En d’autres termes, chaque électron n’interagit pas explicitement avec les autres électrons.

Les trois méthodes HF communément utilisées sont (voir la figure ci-dessous) :

![HF_Orb.png](./Graphics/HF_Orb.png)

1. **RHF (Restricted Hartree – Fock)** utiliśee pour des molécules à couches pleines ou fermées. Les spin-orbitales sont soit $\alpha$, soit $\beta$ et tous les orbitales sont doublement occupés par des spin-orbitales $\alpha$ et $\beta$;</br> 
2. **ROHF (Restricted Open-Shell Hartree–Fock)** utilisée pour des molécules à couches ouvertes où le nombre d'électrons les orbitales n'est pas le même. ROHF utilise autant que possible les orbitales doublement occupées et les orbitales une fois occupées par les électrons non-apariés;
3. **UHF (Unrestricted Hartree-Fock)** utilisée pour des molécules à couches ouvertes où le nombre d'électrons les orbitales n'est pas le même. Les orbitales UHF peuvent avoir avoir des spin $\alpha$ ou $\beta$, mais les orbitales $\alpha$ et $\beta$ peuvent avoir des composants spatiales différents.

Les équivalents DFT sont,

4. **RKS (Kohn-Sham restreint)**,
6. **ROKS (Restricted Open-Shell Kohn – Sham)**,
5. **UKS (Kohn-Sham sans restriction)**.

### PySCF

![Pyscf.png](./Graphics/Pyscf.png)

**Python-based Simulations of Chemistry Framework (PySCF)** est un programme de chimie computationnelle *ab initio* implémenté nativement dans le langage de programme Python. Le package vise à fournir une plate-forme simple, légère et efficace pour le développement et le calcul de codes de chimie quantique. Il fournit diverses fonctions pour faire, au niveau non relativiste, 
* la théorie de Hartree-Fock,
* la Second-order Møller–Plesset perturbation theory (MP2),
* la théorie fonctionnelle de la densité (DFT),
* la Multi-configuration self-consistent field (MCSCF),
* la théorie des clusters couplés (CC).

PySCF est utilisé dans la plupart des packages de calculs quantiques pour les calculs HF (première quantification) et diverses étapes des transformations de seconde quantification.

Il est à noter que la méthode SCF ou Self Consistent Field est une méthode variationnelle où la fonction d'essai est un déterminant de Slater ou une combinaison linéaire de ces déterminants. 

```python
try:
    import pyscf
except:
    %pip install pyscf -U
    import pyscf

pyscf.__version__
```

### Information sur la molécule - Création de l'objet molecule avec `pyscf.gto.Mole()`

Lors de la définition d'une molécule, on spécifie les propriétés moléculaires. 
* L’une d’elles est la charge totale du système, c’est-à-dire la différence entre le nombre d’électrons et la charge nucléaire totale (la valeur par défaut est 0). 
* Une autre quantité importante à spécifier est le spin, défini par la différence entre le nombre d'électrons alpha et bêta (la valeur par défaut est 0). Cela permet de modéliser des systèmes radicalaires : un spin de 0 caractérise une molécule à couche restreinte tandis qu'un spin différent de 0 définit un système à couche ouverte.

Considérons le dimère d'eau (deux molécules d'eau à liaison hydrogène) illustré par la figure suivante

![water_dimer.jpg](./Graphics/water_dimer.jpg)



* Importer les modules nécessaires

```python
from pyscf import gto # Gaussian type orbitals
try :
    import py3Dmol
except:
    %pip install py3Dmol
    import py3Dmol
py3Dmol.__version__
```

* Définir la PySCF molecule object 

```python
# Define the xyz coordinates of the molecule
mol_xyz = """
    O  -1.551007  -0.114520   0.000000
    H  -1.934259   0.762503   0.000000
    H  -0.599677   0.040712   0.000000
    O   1.350625   0.111469   0.000000
    H   1.680398  -0.373741  -0.758561
    H   1.680398  -0.373741   0.758561""" # Ref: https://cccbdb.nist.gov/expgeom1x.asp

smi_key = 'Dwater'

# Define the PySCF molecule object
DWat_mol=gto.Mole(
    atom=mol_xyz,
    basis='def2-TZVP',
    charge=0,      # 0 by default
    spin=0,        # 0 by default, defined as (n_up - n_down == nelec_alpha - nelec_beta)
    unit='Angstrom', # Can also be 'Bohr'
)
DWat_mol.build()

```

```python

"""_3D representation with py3Dmol
"""
# Create a py3Dmol.view object
xyz_view = py3Dmol.view(width=250, height=250)

# Add the molecule model to the view
xyz_view.addModel(DWat_mol.tostring(format="xyz"),'xyz')
# xyz_view.addModel(mol_xyz,'xyz')

# Set the style of representation
xyz_view.setStyle({'stick': {}, 'sphere': {'scale': 0.40}})

# Zoom to fit the molecule in the view
xyz_view.zoomTo()

# Show the molecule view
xyz_view.show()

```

 * Impression de quelques propriétés

```python
print(f'Le nombre total d\'électrons est {DWat_mol.nelectron} \
et le nombre total d\'électrons (alpha, béta) est {DWat_mol.nelec}\n')
print(f'Les index (0-Based) du (HOMO,LUMO) sont {DWat_mol.nelectron//2 -1, DWat_mol.nelectron//2}\n')
print(f'Le nombre d\'orbitales atomiques, dans la base {DWat_mol.basis}, est {DWat_mol.nao}\n')
print(f'L\'énergie nucléaire vaut {DWat_mol.energy_nuc()} Hartrees')
```

Il est à noter que les termes HOMO et LUMO sont utilisés en chimie quantique pour désigner les orbitales moléculaires particulières d'une molécule. 

* **HOMO (Highest Occupied Molecular Orbital)** est la plus haute orbitale moléculaire occupée. Il représente l'orbitale où se trouvent les électrons les plus énergétiques et est donc la plus proche de la couche de valence. Les électrons présents dans l'HOMO sont généralement impliqués dans les réactions chimiques, les liaisons avec d'autres atomes ou molécules et les interactions électroniques.

* **LUMO (Lowest Unoccupied Molecular Orbital)** est la plus basse orbitale moléculaire non occupée. Elle représente l'orbitale qui peut accepter des électrons lors de réactions chimiques ou d'interactions avec d'autres molécules. Le LUMO est généralement d'une énergie plus élevée que l'HOMO et est impliqué dans les processus de transfert d'électrons, de formation de liaisons chimiques ou de réactions de dissociation.

* La différence d'énergie entre l'HOMO et le LUMO, appelée le **gap HOMO-LUMO**, est également un facteur important dans la détermination des propriétés d'absorption et d'émission de lumière des molécules organiques, ce qui est important en chimie des matériaux et en spectroscopie.

<!-- #region -->

### [Hartree-Fock](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method)

* Hatree-Fock (HF) est le point de départ de l'essentiel de la chimie quantique.
* Les orbitales sont optimisés de manière variationnelle pour un seul déterminant de Slater.
* En travaillant sur la base de la fonction de base centrée sur l'atome, on résoud les équations de [Roothaan-Hall](https://en.wikipedia.org/wiki/Roothaan_equations).

$$\mathbf{FC} = \mathbf{SC} \epsilon.$$

* $\mathbf{F}$ est la matrice Fock;
* $\mathbf{C}$ est la matrice des coefficients orbitaux moléculaires;
* $\mathbf{S}$ est la matrice de chevauchement orbitale atomique;
* $\epsilon$ est le vecteur des énergies orbitales moléculaires.

![HF Flowchart](./Graphics/Hartree-Fock.png)

### Création d'un objet mean-field (mf) avec `scf.RHF()`

Pour la plupart des usages, les types de calculs SCF pertinents sont :

1. RHF (Restricted Hartree – Fock)
2. UHF (Hartree-Fock sans restriction)
3. ROHF (Restricted Open-Shell Hartree–Fock)


L'initialisation de l'objet champ moyen pour ces six types de calculs est illustrée ci-dessous. Ces classes nécessitent au moins un argument, à savoir un objet molécule, `mol`.

1. `mf=scf.RHF(mol)`
2. `mf=scf.UHF(mol)`
3. `mf=scf.ROHF(mol)` (où mol est définie comme étant à couche ouverte (open-shell), mol.spin!=0)
4. `mf=scf.RKS(mol)`
5. `mf=scf.UKS(mol)`
6. `mf=scf.ROKS(mol)` (où mol est définie comme open-shell, mol.spin!=0)

* La première étape consiste à créer l'objet mean-field (RHF dans notre cas)
<!-- #endregion -->

```python
from pyscf import scf

mf=scf.RHF(DWat_mol)
```

Par la suite, on exécute le caclul avec la méthode `.run()`ou `.kernel()`.

```python
mf.kernel()
```

Une fois l'exécution du noyau (kernel) SCF achevée, l'objet mean-field (mf) est mis à jour avec plusieurs attributs de sortie utiles :
1. `mf.mo_occ` - Occupation MO (vecteur de longueur égale au nombre de MO)

```python
mf.mo_occ
```

```python
mf.get_occ()
```

```python
try:
    import plotly.express as px
except:
    !pip install plotly
    import plotly.express as px

# Plot the MO Occupations
fig = px.line(y=mf.mo_occ, markers=True, title="Molecular Orbital (MO) Occupations")
fig.update_layout(xaxis_title="Orbital Index (0-Based)", yaxis_title="MO Occupation")
fig.show()
```

On note, sachant que l'indexation commence à 0, que 
* l'orbite moléculaire la plus occupée (HOMO) a pour index 9 (nombres électrons//2 - 1),
* et de l'orbite moléculaire la plus basse inoccupée (LUMO) a pour indexe 10 (nombres électrons//2).


Pour avoir les index du LUMO et l'HOMO à partir du mf:

```python
lumo_idx = mf.mo_occ.tolist().index(0.)
homo_idx = lumo_idx - 1
print(f'Les index (0-Based) du (HOMO,LUMO) sont {homo_idx,lumo_idx}')
```

2. Energies Homo, Lumo et gap en eV

$$E_g = E_{\rm LUMO} - E_{\rm HOMO}$$

```python
from pyscf.data.nist import HARTREE2EV as au2ev 

print(f'Energie du Homo = {mf.mo_energy[homo_idx] * au2ev} eV')
print(f'Energie de Lumo = {mf.mo_energy[lumo_idx] * au2ev} eV')
print(f'Energie du gap Homo-Lumo = {abs(mf.mo_energy[homo_idx] - mf.mo_energy[lumo_idx]) * au2ev} eV')
```

```python
# Plot the MO Energies (i.e. eigenvalues of the Fock matrix)
fig = px.line(y=mf.mo_energy, markers=True, title="Molecular Orbital (MO) Energies (a.u.)")
fig.update_layout(xaxis_title="Orbital Index (0-Based)", yaxis_title="MO Energies (a.u.)")
fig.show()
```

2. `mf.e_tot` - Énergie SCF totale en unités de Hartrees

```python
mf.e_tot
```

2. `mf.dump_scf_summary()` - Synthèse des énergies du calcul SCF effectué.

```python
mf.dump_scf_summary()
```

```python
mf.scf_summary
```

On peut aussi spécifier les paramètres de calculs du mf.

```python
mf=scf.RHF(DWat_mol)

mf.conv_tol=1e-12 # the difference in the SCF energy (in Hartrees) between two successive cycles
mf.conv_tol_grad=1e-8 # the root-mean-square of the orbital gradient
mf.direct_scf_tol=1e-13 #the threshold for discarding integrals
mf.init_guess='atom' #  vital for the efficient completion of the SCF procedure
mf.max_cycle=100 # the maximum number of SCF cycles that should be performed before the calculation terminates. For systems that are notoriously difficult to converge, the value should be increased to 100 or even 1000
mf.max_memory=4000 #determines the maximum amount of memory (in Megabytes) that PySCF is allowed to utilize during the SCF procedure
mf.verbose=5 #controls the print level for the mean-field object (in the 'Dwat_HF.txt' defined above)

mf.kernel()
```

#### Quelques limites des méthodes HF

La méthode de Hartree-Fock présente certaines limites :

 1. **Approximation de l'échange-correlation**. La méthode de Hartree-Fock ne tient pas compte de manière exacte de l'interaction électron-électron appelée **correlation électronique**. Elle utilise une approximation connue sous le nom d'approximation de l'échange-correlation. Cette approximation simplifie l'interaction électronique en la traitant de manière moyennée, ce qui peut conduire à des erreurs significatives dans certains cas.

2. **Non-inclusion de la corrélation électronique dynamique**. La méthode de Hartree-Fock considère les électrons comme un ensemble d'orbitales indépendantes dans un champ moyen créé par les autres électrons. Cependant, elle ne prend pas en compte la corrélation électronique dynamique, c'est-à-dire les effets de corrélation qui dépendent du mouvement des électrons dans le système. Cela limite sa précision pour les systèmes où la corrélation électronique dynamique est importante, tels que les systèmes fortement corrélés ou les réactions chimiques avec des états de transition.

3. **Problèmes avec les systèmes à interaction forte**. La méthode de Hartree-Fock peut être moins précise pour les systèmes contenant des électrons fortement corrélés, tels que les systèmes avec des liaisons multiples ou des électrons non-appariés. Dans ces cas, d'autres méthodes plus avancées, telles que les méthodes de la théorie de la fonctionnelle de la densité (DFT) ou les méthodes de couplage de cluster, peuvent être nécessaires.

4. **Sensibilité aux bases de fonction**. La précision des résultats de la méthode de Hartree-Fock dépend fortement du choix des bases de fonction utilisées pour décrire les orbitales électroniques. Différentes bases de fonction peuvent donner des résultats légèrement différents, et il n'y a pas de base de fonction unique qui convient à tous les systèmes.


<!-- #region -->
### Théorie de la fonctionnelle de la densité - DFT

Dans la [KS-DFT](https://en.wikipedia.org/wiki/Density_functional_theory), proposée pour la première fois par Kohn et Sham en [1964](https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.B864), la densité électronique d'un système de référence sans interaction est utilisée pour représenter la densité du véritable système en interaction. En effet, Kohn et Sham stipulent qu'il existe
1. une cartographie biunivoque entre la densité électronique et le potentiel externe et
2. un principe variationnel pour la densité électronique.

En conséquence, la formulation informatique de KS-DFT ressemble à celle de la théorie de Hartree-Fock (HF), mais avec un potentiel de Fock effectif différent. Dans KS-DFT, l’énergie électronique totale est définie comme suit :

![dft_choices.jpg](./Graphics/dft_choices.jpg)

* $\mathtt{T}_e= -\sum_i\frac{\hbar^2}{2m_e}\nabla^2_i$ est l'énergie cinétique électronique sans interaction;
* $v_{\rm ext} = \int\rho(\mathbf{r})\times\Big[-\sum_{i,I}\frac{Z_I}{|\mathbf{r}_i-\mathbf{R}_I|}\Big]d\mathbf{r}$ est l'énergie due au potentiel externe ou la répulsion Coulombienne entre les électrons et les noyaux;, 
* $ v_H = +\frac12\iint\frac{\rho(\mathbf{r}_i)\rho(\mathbf{r}_j)}{|\mathbf{r}_i-\mathbf{r}_j|}d\mathbf{r}_id\mathbf{r}_j$ est l'énergie Coulombienne; 
* et $v_{\rm XC}$ est l'énergie d'échange-corrélation (*XC*) qui fournit les corrections des termes $\mathtt{T}_e$ et $v_H$. En pratique, $v_{\rm XC}$ est approché par une approximation fonctionnelle de la densité, qui elle-même peuvent être divisée en plusieurs classes le long de différents échelons de l'échelle de Jacob :
    * approximations de densité locale (par exemple LDA ; *XC* l'énergie dépend uniquement de la densité électronique, $\rho$),
    * approximations généralisées du gradient (GGA ; l'énergie *XC* dépend également du gradient de densité, $|\nabla\rho|$),
    * méta-GGA (l'énergie *XC* dépend aussi de la densité d'énergie cinétique et/ou de la densité laplacienne, $\sum_i |\nabla \psi_i|^2$, $\nabla^2\rho$ ; ce dernier n'est pas supporté dans PySCF pour le moment),
    * fonctionnelles de corrélation non locales (l'énergie *XC* implique une double intégrale),
    * fonctionnelles de densité hybride (une fraction de l'échange exact est utilisée), et
    * fonctionnelles de densité corrigées à longue portée (l'échange exact est utilisé avec un noyau d'interaction modifié).
    
![DFT_Accuracy.jpg](./Graphics/DFT_Accuracy.jpg)

Ainsi, l'Hamiltonien à un électron permettant de calculer l'énergie totale d'une orbitale électronique est

$$\hat{h}_{\rm KS} \phi_i(\mathbf{r}_i) = \Big[-\frac12 \nabla_i^2 + -\sum_I\frac{Z_I}{|\mathbf{r}_i-\mathbf{R}_I|} + 
\int\frac{\rho(\mathbf{r}_j)}{|\mathbf{r}_i-\mathbf{r}_j|}d\mathbf{r}_j +  \mathtt{V}_{\rm XC}[\rho(\mathbf{r}_i)] \Big] \phi_i(r) = \epsilon_i \phi_i(\mathbf{r}_i).$$

Les solutions des états propres de cette équation peuvent être utilisées pour trouver une densité qui, par sa construction, reproduit la densité exacte du système complet en interaction :

$$\rho(\mathbf{r}_i)=2\sum_i |\phi_i^{\rm KS}(\mathbf{r}_i)|^2.$$

![kohn_sham_cycle.jpg](Graphics/kohn_sham_cycle.jpg)

### Création d'un objet mean-field (mf) avec `scf.RKS()`

* Exécuter un calcul DFT dans PySCF est presque aussi simple que d'exécuter un calcul Hartree-Fock. 

* PySCF donne aux utilisateurs l'accès à un grand nombre de fonctionnalités via les bibliothèques [libXC](https://tddft.org/programs/libXC/) et [XCFun](https://github.com/dftlibs/XCfun).


Afin d'accéder à des fonctions supplémentaires qui ne relèvent pas de la norme SCF/HF, il est nécessaire d'importer le module dft (théorie fonctionnelle de la densité).
<!-- #endregion -->

```python
from pyscf import dft

mfd=dft.KS(DWat_mol)

mfd.xc='wB97X-D' # la fonctionnelle

mfd.kernel()
```

```python
mfd.dump_scf_summary()
```

```python
mfd.scf_summary
```

Rappellons que précédement avec HF, on avait les énergie suivantes:

```python
mf.dump_scf_summary()
```

RHF et RKS (ainsi que UHF et UKS) partagent de nombreux attributs, nous ne couvrirons donc que ceux qui concernent spécifiquement DFT. Ces attributs supplémentaires concernent principalement la fonctionnelle d'échange-corrélation et le maillage utilisé pour son intégration numérique.


### Analyse du champ moyen

On utilise la méthode `analyze` de l'objet de champ moyen pour imprimer les propriétés moléculaires pertinentes. 

```python
mf.analyze(verbose = 4)
```

#### Quelques limites de la méthode DFT

Bien que la DFT présente de nombreux avantages, elle comporte également certaines limitations importantes :

1. **Approximation de l'échange-correlation**. Tout comme la méthode de Hartree-Fock, la DFT repose sur une approximation de l'échange-correlation. Différentes approximations sont disponibles, telles que l'approximation de la fonctionnelle de l'échange-corrélation locale (LDA) ou l'approximation de la fonctionnelle de l'échange-corrélation généralisée (GGA). Ces approximations peuvent être moins précises dans certains cas, conduisant à des résultats moins précis pour certaines propriétés chimiques.

2. **Difficulté à traiter les interactions à longue portée**. La DFT est moins efficace pour traiter les interactions à longue portée, telles que les interactions de dispersion (van der Waals). Les approximations standard de la DFT ne capturent pas correctement ces interactions, ce qui peut conduire à des erreurs significatives lors de la modélisation de systèmes moléculaires ou de matériaux qui dépendent fortement de ces interactions.

3. **Prise en compte incomplète de la corrélation électronique dynamique**. Bien que la DFT prenne en compte certaines contributions de la corrélation électronique, elle n'inclut généralement pas de manière exacte la corrélation électronique dynamique, qui est importante dans les systèmes avec des états de transition, des réactions chimiques et des propriétés d'excitation. Des méthodes plus avancées, telles que la DFT à deux électrons (DFT à deux particules), sont nécessaires pour capturer plus précisément ces effets de corrélation électronique dynamique.

4. **Sensibilité aux bases de fonction**. Comme pour la méthode de Hartree-Fock, la précision des résultats de la DFT dépend du choix des bases de fonction. Différentes approximations de la fonctionnelle de l'échange-corrélation peuvent nécessiter des bases de fonction spécifiques, et il n'y a pas de base de fonction unique qui convient à tous les systèmes.

5. **Problèmes avec les systèmes fortement corrélés**. La DFT peut être moins précise pour les systèmes avec des électrons fortement corrélés, tels que les systèmes avec des liaisons multiples ou des électrons fortement localisés. Dans ces cas, des méthodes plus avancées, telles que la DFT à la dynamique exacte des échanges, sont nécessaires pour obtenir des résultats précis.



## Choix du set de bases

L'équation de Schrödinger à plusieurs électrons est généralement résolue sur une base d'orbitales moléculaires exprimées sous la forme d'une **combinaison linéaire de ses orbitales atomiques (LCAO, Linear combination of atomic orbitals)**, 
$$
     \phi_a(\mathbf{r}) =
     \sum_A^\mathrm{atomes}
     \sum_\alpha C^A_{\alpha a}
     \chi^A_{\alpha}(\mathbf{r} - \mathbf{R}_A)
$$
Les coefficients LCAO, utilisés pour construire des orbitales moléculaires (MO), sont obtenus grâce à la minimisation de l'énergie totale, selon le principe variationnel. Pour plus de commodité, la base est généralement représentée par des fonctions centrées sur les atomes (les positions atomiques) pour les systèmes moléculaires, ou par des fonctions d'état plane en physique du solide.

Les deux classes d'orbitales de base approximatives couramment utilisées sont :
* les **Slater-type orbitals (STOs)** basées sur le déterminant de Slater, construitent avec des fonctions exponentielles décroissantes $e^{-\alpha r}$,
* et les **orbitales cartésiennes de type Gaussien (GTO)**, construitent avec des fonctions Gaussiennes $e^{-\alpha r^2}$. L’utilisation de ces fonctions Gaussiennes est bien plus efficace que la stratégie STO pour calculer les nombreuses intégrales de chevauchement nécessaires à la chimie quantique.

Ces deux types de fonctions de base peut être combiné comme **STO-nG (Slater-type orbital-n Gaussians)**, où n est le nombre de Gaussiennes utilisées pour faire les approximations.

Les ensembles de base sont construits à partir de nombreuses primitives Gaussiennes. La notation de Pople est l'une des notations les plus explicites pour les ensembles de bases :
- le trait d'union (-) sépare la référence aux orbitales centrales et de valence ;
- les chiffres représentent le nombre de primitives Gaussiennes dans une fonction ;
- il y a autant de chiffres qu'il y a de fonctions pour chaque orbitale atomique (zeta $\zeta$) ;
- la lettre G fait référence à l'utilisation de fonctions Gaussiennes.

À titre d'exemple, l'ensemble de base double zêta `3-21G` utilise trois primitives Gaussiennes pour chaque orbitale centrale et un ensemble de deux plus une primitives Gaussiennes pour les orbitales de valence. Plus d'informations, d'exemples et de discussions sur les GTO peuvent être trouvés dans ce [lien](https://chem.libretexts.org/Courses/Pacific_Union_College/Quantum_Chemistry/11%3A_Computational_Quantum_Chemistry/11.02%3A_Gaussian_Basis_Sets).

Généralement, on adopte des ensembles identiques de fonctions de base pour tous les atomes d'un élément donné et les ensembles adoptés pour les différents éléments doivent avoir une précision équilibrée. Une déclaration telle que « *nous avons utilisé la base 6-31G définie dans le calcul* » signifie qu'un ensemble équilibré d'ensembles d'AO a été adopté pour tous les éléments.

Ci-dessous, des exemples d’ensembles de base populaires utilisés en chimie quantique.


```python
basis_sets = [
    "STO-3G",           # Simple zeta, minimal basis
    "3-21G",            # Double zeta with 3 Gaussian primitives
    "6-31G",            # Double zeta with 6 Gaussian primitives
    "6-31G*",           # Double zeta with 6 Gaussian primitives
    "6-31G(d,p)",       # Polarization functions (+ 5 d-orbitals for all atoms except H, +3 p-orbitals for H atoms) added
    "6-311G(d,p)",      # Triple zeta with polarization functions
    "6-311+G(d,p)",     # Triple zeta with polarization functions and diffuse functions
    "cc-pvqz",          # Quadruple zeta
    "cc-pv5z",          # Quintuple zeta
    "def2-SVP",         # Double zeta with polarization functions
    "def2-SVPD",         # Double zeta with polarization functions and diffuse functions
    "def2-TZVP",        # Triple zeta with polarization functions
    "def2-TZVPD",       # Triple zeta with polarization functions and diffuse functions
    "def2-TZVPP",       # Triple zeta with polarization functions and diffuse functions
    "def2-TZVPPD",      # Triple zeta with polarization functions, diffuse functions, and d functions
    "def2-QZVP",        # Quadruple zeta with polarization functions
    "def2-QZVPP",       # Quadruple zeta with polarization functions and diffuse functions
    "def2-QZVPPD",      # Quadruple zeta with polarization functions, diffuse functions, and d functions
    "def2-QZVPD",       # Quadruple zeta with polarization functions and diffuse functions
    "def2-QZVPPD",      # Quadruple zeta with polarization functions, diffuse functions, and d functions
]
```

Because the basis set is an approximation, it is highly desirable to be able to control its accuracy in order to make tradeoffs between the cost of the calculation and the accuracy of the obtained results. Accordingly, modern basis sets typically come in families of varying size: the smallest sets enable quick but qualitative calculations, while the larger sets enable quantitative computations at the cost of more computer time. In contrast to traditional basis sets, modern basis set families allow for a cost-efficient approach to the complete basis set limit, at which point the error in the one-electron basis set no longer affects the calculation. Note that also other types of basis sets than Gaussians may be used for quantum chemistry.

The Karlsruhe def2 family of Gaussian basis sets,227 which are a good all-round choice for general chemistry as they are available for the whole periodic table up to radon (Z = 86). The Karlsruhe def2 sets come in three levels of accuracy. Split-valence (SV) basis sets are the smallest reasonable basis set for general applications. The def2-SVP basis is a SV basis set with polarization (P) functions, and is similar in size to the 6-31G** also known as the 6-31G(d,p) basis set. Like 6-31G**, the def2-SVP set can also be used without polarization functions on hydrogen atoms; this basis is called def2-SV(P), it is smaller than the 6-31G* basis, and it is often useful for quick qualitative/semi-quantitative calculations. For more quantitative calculations, the def2 series also contains a triple-ζ valence polarization set (def2-TZVP) as well as a quadruple-ζ valence polarization set (def2-QZVP), which typically suffice for achieving the complete basis set limit in HF and DFT calculations. Calculations
at post-HF levels of theory, however, require larger basis sets with additional polarization functions; the def2-TZVPP and def2-QZVPP basis sets exist for this purpose. Diffuse functions (D) are necessary for the proper description of anions as well as to model e.g. electric polarizabilities; sets are likewise available at all levels of accuracy (def2-SVPD, def2-TZVPD, def2-TZVPPD, def2-QZVPD, def2-QZVPPD) for this purpose.

<!-- #region -->
La fonction d'état pour une molécule est une approximation lors de l'utilisation d'un ensemble de base finie,
$$
\Phi_{\text{approx}} = \sum_{i=0}^N c_i |i \rangle
$$

En vertu du principe variationnel, en considérant que $\Phi_{\text{approx}}$ est normalisée, l'énergie approximée est supérieure ou égale à la véritable énergie de l'état fondamentale,
$$
\langle \Phi_{\text{approx}} | \hat{H} | \Phi_{\text{approx}} \rangle = \epsilon_{\text{approx}} \geq  \epsilon_{\text{true}}
$$

A mesure que l'ensemble de base s'agrandit, l'approximation résultante  se rapproche de la véritable fonction d'état,
$$
\Phi_{\text{approx}} = \sum_{i=0}^{\infty} c_i |i \rangle \approx \Psi_{\text{true}}
$$

Cependant, augmenter la taille de l'ensemble de base a pour effet d'augmenter le nombre de paramètres à optimiser, étendant ainsi les ressources de calculs requises dans l'espace (espace de stockage) et en temps (nombre d'opérations). Il est souvent judicieux d'équilibrer la taille de l'ensemble de base pour optimiser le temps de résolution par rapport à la précision cible ($|\epsilon_{\text{approx}} -  \epsilon_{\text{true}}|$), puisque l'augmentation de la complexité de l'ensemble de base produit des rendements décroissants en termes de gain de précision.

Une liste complète des ensembles de bases disponibles pris en charge dans [PySCF](https://pyscf.org/) est disponible sur ce [lien](https://pyscf.org/_modules/pyscf/gto/basis.html). De plus en plus de bases sont disponibles pour la communauté, et le choix d'un ensemble de bases peut s'avérer une tâche difficile. La théorie stipule qu'une base infinie est le moyen d'obtenir de vrais calculs des fonctions d'état. En pratique, il n'est pas possible d'effectuer les calculs avec des ensembles de bases infinis et d'obtenir des vraies fonction d'état. On peut cependant sélectionner des ensembles de bases finis et adapter nos approches aux cas d'utilisation. En règle générale, 
* un double zêta (par exemple `6-31G` ou `cc-pvdz`) est la taille de base minimale requise pour obtenir des énergies semi-quantitatives;
* pour les molécules présentant une différence d'électronégativité élevée entre les atomes, l'ajout de fonctions de polarisation (par exemple `6-31G**` ou `6-31G(d,p)`) est recommandé;
* les fonctions diffuses (par exemple `6-31G+(d,p)`) doivent être incluses lorsqu'il s'agit d'anions. Les fonctions diffuses sont des fonctions avec de petits coefficients (petit $\alpha$ dans $e^{-\alpha r^2}$) qui peuvent s'adapter à la densité électronique loin du noyau;
* les ensembles de bases des potentiels de noyau effectif (Effective Core Potential, ECP) sont nécessaires lorsqu'il s'agit d'atomes plus lourds avec des effets relativistes et des effets non triviaux (par exemple, l'atome d'iode est souvent modélisé dans l'ensemble de base `LANDL2DZ`).


**Tableau** : Catégories d'ensemble de base.
| Catégorie | Caractère | Exemple |
| -------- | --------- | ------- |
| Minimale| Fonction de base unique par orbitale atomique occupée, la taille des éléments de la deuxième rangée est $[2s1p]$ | STO-3G |
| Double-$\zeta$ | Deux fonctions de base par orbitale atomique de valence occupée, la taille des éléments de la deuxième rangée est $[3s2p]$ | 6-31G |
| Triple-$\zeta$ | Trois fonctions de base par orbitale atomique de valence occupée, pour les éléments de deuxième rangée $[4s3p]$ | 6-311G |
| Polarisation | Fonctions de base avec un moment cinétique plus élevé ; important pour la liaison chimique ; la taille des éléments de la deuxième ligne est $[3s2p1d]$ | cc-pVDZ |
| Diffus | Fonctions de base avec de petits exposants ; important pour les états excités ; les tailles des éléments de la deuxième rangée sont $[4s3p]$ pour 6-31+G et $[4s3p2d]$ pour aug-cc-pVDZ | 6-31+G, aug-cc-pVDZ |

<!-- #endregion -->

### Effet de l'ensemble de base sur les calculs d'énergie totale et le temps de calculs

```python
import time

mf_energies = list()
mf_times = list()
nb_prim = list()

# Perform a Mean-Field calculation for each basis set
for bs in basis_sets:

    # Measure execution time
    start = time.time()
    DWat_mol.basis = bs
    DWat_mol.build()

    mf = scf.RHF(DWat_mol)
    mf.kernel()
    end = time.time()

    nb_prim.append(DWat_mol.npgto_nr())
    mf_energies.append(mf.e_tot)
    mf_times.append(end-start)

```

```python
# Create the results dataframe
import pandas as pd

df_HF = pd.DataFrame({"Basis":basis_sets, 
                      'Nb of GTO primitives':nb_prim, 
                      'Total energy':mf_energies,
                      "Time":mf_times})

df_HF

```

* Illustration avec `matplotlib.pyplot` qui est par défaut statique, mais que nous allons rendre interactif avec la commande magique `%matplotlib notebook`.

```python
%matplotlib notebook
import matplotlib.pyplot as plt

# Create the matplotlib figure
fig, ax = plt.subplots(figsize=(8,5))

# Plot the energies.
ax.set_xticks(range(len(basis_sets)), basis_sets, rotation=90)
ax.set_xlabel("Basis set")
ax.set_ylabel("Energy(Hartree)", color="b")
ax.scatter(range(len(basis_sets)), mf_energies, marker="o", s=50, color="b")

# Plot the time to solution
ax_time = ax.twinx()
ax_time.scatter(range(len(basis_sets)), mf_times, marker="s", s=50, color="r")
ax_time.set_ylabel("Time to solution(s)", color="r", rotation=270, va="bottom")

# Show the graph
plt.tick_params(axis="both", direction="in")
plt.show()
```

* Illustration avec `plotly.express` qui est interactif par défaut.

```python
import plotly.express as px

# Create the figure
fig = px.scatter()

# Convert basis_sets to a list
basis_indices = list(range(len(basis_sets)))

# Plot the energies
fig.add_scatter(x=basis_indices, y=mf_energies, mode="markers",\
                marker=dict(symbol="circle", size=10, color="blue"), name="Energy (Hartree)")

# Plot the time to solution
fig.add_scatter(x=basis_indices, y=mf_times, mode="markers",\
                marker=dict(symbol="square", size=10, color="red"), name="Time to solution (s)", yaxis="y2")

# Set the x-axis labels and tick labels
fig.update_layout(xaxis=dict(tickmode="array", tickvals=basis_indices, ticktext=basis_sets, tickangle=90))

# Set the y-axis labels and colors
fig.update_layout(yaxis=dict(title="Energy(Hartree)", color="blue"),\
                  yaxis2=dict(title="Time to solution (s)", color="red", overlaying="y",\
                              side="right"))

# Show the graph
fig.show()
```

Le tableau précédent et cette figure suggèrent `def2-TZVPP` comme le bon ensemble de bases pour cette molécule.


## Quelques méthodes post-HF

Dans ce qui suit, pour des raisons computationnelles, nous allons utiliser l'ensemble de base `STO-3G` et la molécule `LiH`.

### Théorie des perturbations de Møller-Plesset - MP2

* La [MP2](https://en.wikipedia.org/wiki/M%C3%B8ller%E2%80%93Plesset_perturbation_theory) effectue des corrections perturbatives de l'approximation Hartree-Fock.

```python
from pyscf import scf, gto
import time

# Mean-field HF
mol = gto.M(
    atom = "Li 0 0 0; H 0 0 1.5474",
    basis = 'sto-3g',
)

start = time.time()
myhf = scf.RHF(mol).run()
myhf_time = time.time() - start
```

```python
from pyscf import mp

start = time.time()
mymp2 = mp.MP2(myhf).run()
mymp2_time = time.time() - start + myhf_time
```

#### Cluster couplé - CC

* La [CC](https://en.wikipedia.org/wiki/Coupled_cluster) est une méthode perturbative qui améliore l'approximation de Hartree-Fock.
* Les clusters couplés simples et doubles (CCSD) incluent une excitation simple et double en plus de la fonction d'état HF.
* La précision peut être améliorée en incluant les triples de manière perturbatrice (CCSD(T)).
* Description non variationnelle, mais détaillée des états fondamentaux. 
* L'extension aux états excités est l'EOM-CCSD.

```python
from pyscf import cc

start = time.time()
mycc = cc.CCSD(myhf).run()
mycc_time = time.time() - start + myhf_time

e_ccsd_t = mycc.ccsd_t()
mycct_time = time.time() - start + myhf_time
```

### Interaction de configuration complète - FCI

* L'[interaction de configuration complète (FCI)](https://en.wikipedia.org/wiki/Full_configuration_interaction) est exacte pour un choix donné d'ensemble de base.
* FCI inclut tous les déterminants de Slater de symétrie appropriée dans la base du vecteur propre.
* Le coût computationnel augmente de façon exponentielle avec la taille du système.
* Également connue sous le nom de **diagonalisation exacte**.

```python
from pyscf import fci

start = time.time()
myfci = fci.FCI(myhf).kernel()
myfci_time = time.time() - start 
```

### Synthèse des méthodes HF et post-HF

Les données importantes sont enregistrées dans les objets de la méthode PYSCF, ce qui facilite l'analyse et la visualisation.

```python
# Collect data

methods = ["HF", "MP2", "CCSD", "CCSD(T)", "FCI"]
energies = [myhf.e_tot, mymp2.e_tot, mycc.e_tot, mycc.e_tot + e_ccsd_t, myfci[0]]
mf_times = [myhf_time, mymp2_time, mycc_time, mycct_time, myfci_time]

# Create the results dataframe
import pandas as pd

df_Eies = pd.DataFrame({"Energy (a.u.)":energies, "Time (s)":mf_times, 
                        "Delta E(kcal/mol)": [(abs(energy - myfci[0]) * 627.509474) for energy in energies]}, 
                       index = methods)

df_Eies
```

On note que les méthodes HF et MP2 produisent ici des différences d'énergie12.3 et 4.5 kcal/mol respectivement, bien supérieure à [précision chimique](https://en.wikipedia.org/wiki/Computational_chemistry) (1 kcal/mol).

```python
import plotly.express as px

# Create the figure
fig = px.scatter()

# Plot the energies
fig.add_scatter(x=methods, y=energies, mode="markers",\
                marker=dict(symbol="circle", size=10, color="blue"), name="Energy (a.u.)")

# Plot the time to solution
fig.add_scatter(x=methods, y= mf_times, mode="markers",\
                marker=dict(symbol="square", size=10, color="red"), name="Time to solution (s)", yaxis="y2")

# Set the x-axis labels and tick labels
fig.update_layout(xaxis=dict(tickmode="array", tickvals=methods, ticktext=methods, tickangle=90))

# Set the y-axis labels and colors
fig.update_layout(yaxis=dict(title="Energy(a.u.)", color="blue"),\
                  yaxis2=dict(title="Time to solution (s)", color="red", overlaying="y",\
                              side="right"))

# Show the graph
fig.show()
```

De ce qui précède, il apparaît que, pour un même ensemble de bases chimiques de dimension moyenne, de taille $N$,
* en termes de **précision**, on a l’ordre suivant :

HF < MP2 < CCSD < CCSD(T) < FCI

* tn termes de **temps de calculs**, 

HF ($N^4$) < MP2 ($N^5$) < CCSD ($N^6$)  < CCSD(T) ($N^7$) < FCI ($N!$) 
On note sur ce graphique que l’augmentation de la complexité de l'ensemble de base rapproche la valeur propre de l’énergie de la véritable énergie de Hartree-Fock. Cette tendance peut également être généralisée à d'autres méthodes, comme [CCSD](https://en.wikipedia.org/wiki/Coupled_cluster), [MP2](https://en.wikipedia.org/wiki/M%C3%B8ller%E2%80%93Plesset_perturbation_theory), [FCI](https://en.wikipedia.org/wiki/Full_configuration_interaction), etc., bien que le temps de résolution serait plus long pour ces approches post Hartree-Fock. Ce qu’il faut retenir, c’est que même si une base infinie est nécessaire pour converger vers la véritable énergie, il y a des rendements décroissants à un moment donné, surtout si l’on considère le temps requis pour le calcul. Dans notre cas, `6-311G(d,p)` ou `6-311+G(d,p)` apparaît comme un bon choix en termes de précision et de temps de calcul.

![HF_Post_HF](./Graphics/HF_Post_HF.png)

### Prise en compte de l'énergie de corrélation - Potentiel de l'informatique quantique

Les méthodes Post Hatree-Fock sont destinées à prendre en compte l'énergie de corrélation électronique. L'énergie de corrélation est définie comme
$$
E_{\text{corr}} = E_{\text{FCI}} - E_{\text{HF}},
$$
dans lequel $E_{\text{FCI}}$ et $E_{\text{HF}}$ sont les énergies calculées à partir de la diagonalisation exacte de l'Hamiltonien moléculaire dans l'approximation de Born-Oppenheimer, c'est-à-dire l'interaction de configuration complète (FCI), et la méthode Hartree-Fock, respectivement. Puisqu'il est prohibitif sur le plan computationnel de résoudre le problème FCI, à l'exception de très petits systèmes, il est courant d'introduire des approximations pour inclure les effets de corrélation dans les calculs de structure électronique. Le problème FCI considère la diagonalisation de l'Hamiltonien dans une base de tous les déterminants excités possibles, c'est-à-dire tous les déterminants possibles à 1, 2, 3,…, N corps pour un système à N électrons. Par conséquent, la fonction d’onde FCI est une combinaison linéaire de tous les déterminants possibles pour un système à N électrons. Les méthodes approximatives n'incluent qu'un sous-ensemble des déterminants excités dans la fonction d'état approximative afin de capturer les effets de corrélation les plus importants (par exemple, théorie des perturbations à plusieurs corps, théorie des clusters couplés).

Cette énergie de corrélation peut être divisée en deux catégories : **corrélation statique** et **dynamique**. Pour considérer le premier cas, il faut utiliser des approches multi-configurations. Ces dernières peuvent être abordées avec des méthodes post Hartree-Fock avec un seul état de référence, comme Configuration Interaction (CI), Coupled-Cluster (CC) ou des méthodes perturbatives (MP2, MP3, .etc).
