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

<!-- #region id="yqbzsX3L6AH5" -->
# PHY4268 Tutorial 4 - Qiskit VQE with active space

1. **S. G. Nana Engo**, serge.nana-engo@facsciences-uy1.cm
    * Department of Physics, Faculty of Science, University of Yaounde I
2. **J-P. Tchapet Njafa**, jean-pierre.tchapet-njafa@univ-maroua.cm
    * Department of Physics, Faculty of Science, University of Maroua
       
May 2024
<!-- #endregion -->

<!-- #region id="gph6lvv66AH8" -->
A la fin de ce tutoriel, l'apprenant doit être capable d'utiliser Qiskit-nature pour
1. calculer l'énergie de l'état fondamental d'un Hamiltonien quelconque;
2. définir l'espace actif d'un système moléculaire autour du niveau de Fermi (Niveau HOMO-LUMO);
3. calculer les énergies de l'état fondamental et des premiers états excités y afférent.

<!-- #endregion -->

<!-- #region id="X1UkxnGt6AH9" -->
## Hamiltonien quelconque

On considère l'Hamiltonien
$$\mathtt{H} = 0.4\mathbb{I}\mathtt{X} + 0.6\mathbb{I}\mathtt{Z} + 0.8\mathtt{XY}$$

Pour un $|\psi\rangle$ donné nous voulons évaluer la valeur moyenne de cet Hamiltonien :

$$\langle \mathtt{H} \rangle = \langle \psi |\mathtt{H}| \psi \rangle = 0.4 \langle \psi | \mathbb{I}\mathtt{X} |\psi \rangle + 0.6 \langle \psi | \mathbb{I}\mathtt{Z} | \psi \rangle + 0.8 \langle \psi | \mathtt{XY} | \psi \rangle.$$

Comme on peut voir la valeur moyenne $\langle\mathtt{H} \rangle$ pourrait être calculée en ajoutant les valeurs moyennes de ses parties (termes de Pauli). L'algorithme fait exactement cela. Il construit un circuit quantique pour chaque terme de Pauli et calcule la valeur moyenne du terme de Pauli correspondant. Ensuite, l'algorithme additionne toutes les valeurs moyennes calculées des termes de Pauli et obtient la valeur moyenne de $\mathtt{H}$.

<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 35} id="knZ2hyJg6R9b" outputId="7a1ddb8d-034c-4128-b4a5-2641de50b551"
try:
    import qiskit_nature
except:
    %pip install qiskit-nature[pyscf] -U
    import qiskit_nature

qiskit_nature.__version__
```

```python colab={"base_uri": "https://localhost:8080/"} id="kT1jbfA36AH-" outputId="ebb027b1-f182-4eff-bcfb-cd30f60ad92b"
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import NumPyMinimumEigensolver

H = SparsePauliOp.from_list([("IX", 0.4), ("IZ", 0.6), ("XY", 0.8)])
print(H)

exact_result = NumPyMinimumEigensolver().compute_minimum_eigenvalue(H)
print(f'\n The exact ground state energy is: {exact_result.eigenvalue}')
```

<!-- #region id="2Ookb3_W6AH_" -->
Si on a un Hamiltonien
$$\mathtt{H} = a\mathbb{I} + b\mathtt{Z} + c\mathtt{X} + d\mathtt{Y},\qquad a,b,c,d\in\mathbb{R},$$
il faut faire appel à des nombres arbitraires pour a, b, c et d.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="HNgzJPF56AH_" outputId="f9daed0b-5546-4247-b3ba-a1d27cc00207"
from numpy import random

a, b, c, d = 10*random.sample(4)

H = SparsePauliOp.from_list([("I", a), ("Z", b), ("X", c), ("Y", d)])

print(H)
exact_result = NumPyMinimumEigensolver().compute_minimum_eigenvalue(H)
print(f'\n The exact ground state energy is: {exact_result.eigenvalue}')
```

<!-- #region id="Zz_zVmc26AIA" -->
## Molécule d'eau

### Hamiltionien électronique

![Water_HF_references](Graphics/Water_HF_references.png)

La molécule d'eau a un nombre total de dix électrons. Le déterminant de Slater résultant de l'occupation des cinq premières orbitales moléculaires d'énergie la plus basse avec deux électrons *appariés* dans chaque orbitale, l'un avec spin-up et l'autre avec spin-down, est dit être un état HF à couche fermée avec une *multiplicité de spin=1*. Alternativement, si nous définissons une occupation où les quatre premières orbitales sont doublement occupées et les deux suivantes sont occupées individuellement par des électrons *non appariés* avec spin-up, on dit qu'il s'agit d'un état HF à couche ouverte avec une *multiplicité de spin=3*.

Il n'est pas superflue de rappeler à nouveau ce qui suit.

* La multiplicité, que nous pouvons définir comme $(N_{\rm unpaired}^e + 1)$ avec $N_{\rm unpaired}^e$ le nombre d'électrons non appariés, détermine l'occupation des orbitales moléculaires dans les calculs HF.

* Les orbitales moléculaires sont généralement représentées comme une combinaison linéaire de
**orbitales atomiques**. Les coefficients de d'expansion dans la base atomique sont
calculé en utilisant la méthode de Hartree-Fock (HF). Dans l'approximation HF, chaque électron de la molécule est traité comme une particule **indépendante** qui se déplace sous l'influence du Coulomb potentiel dû aux noyaux, et un champ moyen généré par tous les autres
électrons. Les coefficients optimisés sont précisément ce dont on a besoin pour
construire le Hamiltonien de la seconde quantification.

### Active-Space reduction

En général, les méthodes de corrélation d'électrons post-Hartree-Fock étendent la fonction d'état de la molécule autour de la solution Hartree-Fock, en ajoutant des déterminants de Slater, communément appelés **configurations**, qui résultent de l'excitation des électrons des orbitales HF occupées vers les orbitales inoccupées. Malgré le fait qu'il existe différentes techniques pour tronquer cette expansion, le nombre de configurations augmente de manière combinatoire avec le nombre d'électrons et de fonctions de base et la tâche de trouver les coefficients d'expansion de la fonction d'état devient numériquement insoluble si nous voulons inclure l'ensemble complet des orbitales moléculaires. .

Afin de contourner l'explosion combinatoire, nous pouvons créer un espace actif en classant les orbitales moléculaires en orbitales doublement occupées, actives et externes :

* les orbitales doublement occupées (**core orbitals**) sont toujours occupées par deux électrons;
* les orbitales actives (**valence orbitals**) peuvent être occupées par zéro, un ou deux électrons.
* les orbitales externes (**virtual orbitals**) ne sont jamais occupées.

![Sketch_active_space](Graphics/Sketch_active_space.png)

Dans cette approximation, un certain nombre d'*électrons actifs* peuvent peupler les *orbitales actives* à partir desquelles nous pouvons générer un espace de taille finie de déterminants de Slater.

**Note**</br>
Le nombre de *spin-orbitales actives* détermine le *nombre de qubits* requis pour effectuer des simulations quantiques de la structure électronique de la molécule.

Pour le cas de la molécule d'eau décrite à l'aide d'un ensemble de base minimal, nous avons un total de dix électrons occupant les cinq premières des sept orbitales moléculaires dans l'état de référence HF.

Avec `qiskit_nature.second_q.transformers.ActiveSpaceTransformer`, la réduction se fait en calculant l'opérateur de Fock inactif qui est défini comme
$$ F^I_{pq} = h_{pq} + \sum_i (2 g_{iipq} - g_{iqpi}),
$$
et l'énergie inactive qui est donnée par
$$ E^I = \sum_j h_{ji} + F^I_{jj} = \frac12 \sum_{ij} \Big(h_{ij} + F^I_{jj}\Big) D^I_{ij} ,$$
où $i$ et $j$ itèrent sur les orbitales inactives. En utilisant l'opérateur de Fock inactif à la place des intégrales à un électron, la description de l'espace actif contient un potentiel effectif généré par les électrons inactifs. Par conséquent, cette méthode permet l'exclusion des électrons non centraux tout en conservant une description de haute qualité du système.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 35} id="jYpuxgfi6AIB" outputId="302e9176-6c69-4714-eeb7-60de5ca01616"
try:
    import qiskit_nature
except:
    %pip install qiskit-nature[pyscf] -U
    import qiskit_nature

qiskit_nature.__version__
```

```python id="PUNh3u3P6AIB"
from qiskit_nature.second_q.drivers import PySCFDriver

H2O_driver = PySCFDriver(
    atom='O 0.0 0.0 0.0; H 0.757 0.586 0.0; H -0.757 0.586 0.0',
    basis="sto-3g"
)

H2O_problem = H2O_driver.run()
```

```python colab={"base_uri": "https://localhost:8080/"} id="2qITYn7B6AIC" outputId="a7f997f8-b81d-40c1-b12e-444bd92d4bf1"
print(f'Molecule strcuture: {H2O_problem.molecule}')
print(f'Reference energy: {H2O_problem.reference_energy}')
print(f'Nuclear repulsion energy: {H2O_problem.nuclear_repulsion_energy}')
print(f'Number of particules: {H2O_problem.num_particles}')
print(f'Number of spatial orbitals: {H2O_problem.num_spatial_orbitals}')
print(f'Number of molecular orbitals: {H2O_problem.num_spin_orbitals}')
print(f'Number of alpha electrons: {H2O_problem.num_alpha}')
print(f'Number of beta electrons: {H2O_problem.num_beta}')
print(f'Occupations of the alpha-spin orbitals: {H2O_problem.orbital_occupations}')

```

```python colab={"base_uri": "https://localhost:8080/"} id="zBT-ya926AID" outputId="984e29e5-2acc-4a95-e60e-7ba30097cfed"

print(f"HOMO and LUMO indexes are (Fermi level) {H2O_problem.num_alpha-1, H2O_problem.num_alpha}")
print(f"HOMO and LUMO energies are {H2O_problem.orbital_energies[H2O_problem.num_alpha-1],H2O_problem.orbital_energies[H2O_problem.num_alpha]}")
print(f"HOMO-LUMO gap is {abs(H2O_problem.orbital_energies[H2O_problem.num_alpha-1]-H2O_problem.orbital_energies[H2O_problem.num_alpha])}")

```

<!-- #region id="kPV0y6Rr6AID" -->
### 1.9.3. <a id='toc1_9_3_'></a>[`FreezeCoreTransformer`](#toc0_)

Ce transformateur vous offre un moyen très simple de geler les _core orbitales_ du système moléculaire. Il nécessite que votre problème contienne l'attribut `.molecule` à partir duquel il peut extraire les informations atomiques nécessaires pour effectuer cette réduction d'espace de Hilbert.

Appliquons `qiskit_nature.second_q.transformers.FreezeCoreTransforme` à notre molécule qui, dans ce cas, supprimera la seule orbitale d'énergie la plus basse (réduisant le nombre total d'orbitales spatiales de 7 à 6) et supprimant également les deux électrons de l'intérieur de cette orbitale (comme reflété par le nombre modifié de particules).
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="XbMuEZ9-6AIE" outputId="79263917-5ae2-4d46-e46e-fcfcef2c3d8a"
from qiskit_nature.second_q.transformers import FreezeCoreTransformer

fc_transformer = FreezeCoreTransformer()

fc_H2O_problem = fc_transformer.transform(H2O_problem)
print(f'Number of particules with FC: {fc_H2O_problem.num_particles}')
print(f'Number of spatial orbitals with FC: {fc_H2O_problem.num_spatial_orbitals}')

```

<!-- #region id="8rVsAyQw6AIE" -->
Notez que cette transformation se traduira par un décalage d'énergie constant résultant de la suppression des électrons du noyau. Ce décalage est enregistré à l'intérieur de l'attribut `constants`" de l'Hamiltonien.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="rFbfFtHR6AIE" outputId="feebc477-3882-48b5-c6c8-be3a59efa735"
print(fc_H2O_problem.hamiltonian.constants)
```

<!-- #region id="OnIAlGVM6AIE" -->
Vous pouvez fournir une liste d'indices orbitaux (se rappeler qu'on commence par 0) qui doivent être supprimés du système.

> **Remarque :** Ces orbitales *doivent* être inoccupées, sinon vous subirez une erreur importante dans votre calcul. Vous devez absolument savoir quelles orbitales vous supprimez, car la suppression des mauvaises orbitales peut toujours entraîner de grandes erreurs si la dynamique des systèmes est modifiée de manière significative.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="O8N7fnnU6AIF" outputId="dd1f84f1-fb5b-4a8b-879b-66ca78fc6bfa"
fc_H2O_transformer = FreezeCoreTransformer(remove_orbitals=[5, 6])

fc_H2O_problem = fc_H2O_transformer.transform(H2O_problem)
print(fc_H2O_problem.num_particles)
print(fc_H2O_problem.num_spatial_orbitals)
```

<!-- #region id="QajoZ7lr6AIF" -->
### `ActiveSpaceTransformer`

Ce transformateur généralise la réduction de l'espace de Hilbert qui est effectuée par le `FreezeCoreTransformer`. En termes simples, il vous permet de spécifier un _active space_ en sélectionnant le nombre d'électrons actifs et le nombre d'orbitales spatiales actives. Selon ces paramètres, l'espace actif sera choisi autour du niveau de Fermi.

L'espace actif peut être configuré de l'une des manières suivantes via l'initialiseur :

- lorsque seuls ``num_electrons`` et ``num_spatial_orbitals`` sont spécifiés, ces entiers indiquent respectivement le nombre d’électrons actifs et d’orbitales. L'espace actif sera puis être choisi autour du niveau de Fermi, ce qui donne un choix unique pour n'importe quelle paire de nombres. Néanmoins, les critères suivants doivent être remplis :

   * le nombre restant d'électrons inactifs doit être un nombre pair et positif;
   * le nombre d'orbitales actives ``num_spatial_orbitals`` ne doit pas dépasser le nombre total d'orbitales moins le nombre d'orbitales occupées par les électrons inactifs;

- lorsque ``num_electrons`` est un tuple, cela doit indiquer le nombre de spins alpha et bêta électrons, respectivement. Les mêmes exigences que celles énumérées précédemment doivent être remplies;
- enfin, il est possible de sélectionner un ensemble personnalisé d'orbitales actives via leurs indices en utilisant ``active_orbitals``. Ceci permet de sélectionner un espace actif qui n'est pas placé autour du niveau de Fermi comme décrit dans le premier cas ci-dessus. Lorsque vous utilisez cet argument de mot-clé, les critères suivants doivent être remplis *en plus* de ceux énumérés ci-dessus :

   * la longueur de `active_orbitals` doit être égale à ``num_spatial_orbitals``. Noter que qiskit-nature de  **déduit pas** le nombre d'orbitales actives à partir de cette liste d'indices !

   * lors de l'utilisation d'un tuple de listes pour indiquer les indices orbitaux de spin alpha et bêta séparément, les deux listes doivent remplir le critère précédent;

   * le plus grand indice orbital ne peut **pas** dépasser le ``num_spatial_orbitals`` disponible.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="4vauSHPr6AIG" outputId="9cad8d3a-3136-49c5-8883-79c4a454676b"
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

# We want to reduce it to an active space of 4 electrons in 4 orbitals (HOMO-1,HOMO, LUMO, LUMO+1)
as_transformer = ActiveSpaceTransformer(4, 4) # ActiveSpaceTransformer(num_electrons, num_spatial_orbitals)

H2O_as_problem = as_transformer.transform(H2O_problem)
print(H2O_as_problem.num_particles)
print(H2O_as_problem.num_spatial_orbitals)
```

<!-- #region id="VSbMQyLz6AIG" -->
`qiskit_nature.second_q.transformers.ActiveSpaceTransformer` permet aussi de spécifier manuellement les indices des orbitales actives. Cela vous permet de sélectionner manuellement des _active space_ qui ne se trouvent pas en permanence autour du niveau de Fermi.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="dA4HVZ1z6AIG" outputId="b8bc1e7a-3532-46c1-ea54-637493715d7c"
as_transformer = ActiveSpaceTransformer(2, 2, active_orbitals=[3, 6])

H2O_as_problem = as_transformer.transform(H2O_problem)
print(H2O_as_problem.num_particles)
print(H2O_as_problem.num_spatial_orbitals)
```

<!-- #region id="UGyaGjhs6AIG" -->
### Energie fondamentale de la molecule 1,2-Thiazole atom

Comme vous le savez déjà, il est difficile d'obtenir le Hamiltonien qubit de cette molecule, au vue de sa taille. Nous allons donc utiliser `ActiveSpaceTransformer` pour réduire sa taille.
<!-- #endregion -->

```python id="r03wCpcW6AIG"
# Molecular structure
driver = PySCFDriver(
    atom="""C      1.1291      0.0795     -0.5259
  C      0.7115     -1.2207     -0.4748
  H      2.0789      0.4171     -0.9381
  H      1.2719     -2.0822     -0.8377
  S     -0.0500      1.1306      0.1514
  N     -1.1147     -0.1822      0.5074
  C     -0.5926     -1.3347      0.1299
  H     -1.1422     -2.2662      0.2851""", #C1,2-Thiazole atom
    basis='STO-3G'
)
# Electronic structure problem
problem = driver.run()
```

```python colab={"base_uri": "https://localhost:8080/"} id="E8yK4vgU6AIH" outputId="636a0456-4e74-4da9-f69d-9ca34fd775c8"
# Some properties
print(f"Molecule, basis: {driver.basis}, Hartree-Fock calculation")
print(f"Number of alpha electrons: {problem.num_alpha}")
print(f"Number of beta electrons: {problem.num_beta}")
print(f"Number of spin orbitals: {problem.num_spin_orbitals}")
print(f"Spin orbitals occupation: {problem.orbital_occupations}")
print(f"Spin orbitals energies: {problem.orbital_energies}")
print(f"Molecule reference (HF) total energy: {problem.reference_energy} Ha")
print(f"Molecule nuclear repulsion energy: {problem.nuclear_repulsion_energy} Ha")

```

```python colab={"base_uri": "https://localhost:8080/"} id="JsZDZxHj6AIH" outputId="545ebe0b-8d76-4980-f373-0d61a107b3f3"
# Fermi Level and properties
from pyscf.data import nist
au2ev = nist.HARTREE2EV

print(f"HOMO and LUMO index are {problem.num_alpha-1, problem.num_alpha}")
print(f"HOMO and LUMO energies in eV are\
{problem.orbital_energies[problem.num_alpha-1],problem.orbital_energies[problem.num_alpha] * au2ev}")
print(f"HOMO-LUMO gap in eV is \
{abs(problem.orbital_energies[problem.num_alpha-1]-problem.orbital_energies[problem.num_alpha]) * au2ev}")
```

```python colab={"base_uri": "https://localhost:8080/"} id="a_TUi2ER6AIH" outputId="ac0a82c9-64be-4d2e-ec58-6a1fee6f2e4a"
problem.num_particles
```

<!-- #region id="Fo8hh3RJ6AIH" -->
* Définition d'un espace actif autour du niveau de Fermi

Nous allons choisir un espace actif avec 4 orbitale spatiale, c'est-à-dire (HOMO-1,HOMO, LUMO, LUMO+1).
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="jPukFijq6AII" outputId="81bb0948-de42-471b-e494-874b7f3c6ceb"
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
problem = ActiveSpaceTransformer(4, 4).transform(problem) # Utiliser (2,2) pour rendre le calcul moins couteux (+perte précision)problem.num_particles
problem.num_particles
```

<!-- #region id="qY44qdww6AII" -->
On note que nous sommes partis de (22,22) à (2,2) particules alpha et beta!
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="W1oCXNNo6AII" outputId="5cd97039-7efe-4e23-df1d-3a31b032169f"
hamiltonian = problem.hamiltonian # Hamiltonien en 1ere quantification
fermionic_op = hamiltonian.second_q_op() # Hamiltonien fermionique
if len(fermionic_op) <= 20:
    print(fermionic_op)
else: # print the first 20 terms of the fermionic Hamiltonian operator of the molecule
    print("\n".join(str(fermionic_op).splitlines()[:22] + ["..."]))
```

<!-- #region id="UwKNfwyH6AII" -->
* Hamiltonien qubit avec la réduction $\mathbb{Z}_2$
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="yon9s9s56AII" outputId="b11201d0-eccd-4162-f77f-e7bdfeb2c15a"
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper, BravyiKitaevMapper, TaperedQubitMapper

mapper = ParityMapper(num_particles=problem.num_particles)
mapper = problem.get_tapered_mapper(mapper)
Hamil_z2qubit = mapper.map(fermionic_op)

print(f"Number of items in the PM Z2 Pauli list:", len(Hamil_z2qubit))
if len(Hamil_z2qubit) <= 10:
    print(Hamil_z2qubit)
else:
    print(Hamil_z2qubit[0:10])
```

<!-- #region id="4vUYVu_p6AIJ" -->
* Circuit de l'état initial
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 217} id="08E7-yxL6AIK" outputId="ab10ce03-ab24-451a-f895-f45d728ad8dc"
from qiskit_nature.second_q.circuit.library import HartreeFock

HF_state = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)
HF_state.draw('text',initial_state=True)
```

<!-- #region id="lab12QkT6AIK" -->
Ce qui précède est une solution au **Devoir sur le Hamiltonien Qubit**.

* Définissons le solveur VQE.
<!-- #endregion -->

```python id="KjDQJ93f6AIK"
from qiskit_nature.second_q.circuit.library import UCCSD
from qiskit.primitives import Estimator
from qiskit_algorithms.optimizers import SLSQP
from qiskit_algorithms import VQE
import numpy as np

ansatz = UCCSD(
    problem.num_spatial_orbitals,
    problem.num_particles,
    mapper,
    initial_state = HF_state)

vqe_solver = VQE(Estimator(), ansatz, SLSQP())
vqe_solver.initial_point = np.zeros(ansatz.num_parameters) # Initialisation with zero-qubit
```

<!-- #region id="_KTit7UN6AIK" -->
* Calcul et résultats
<!-- #endregion -->

```python id="rc8WzyO26AIL"
from qiskit_nature.second_q.algorithms import GroundStateEigensolver

# Résultats VQE
GS_VQE = GroundStateEigensolver(mapper, vqe_solver)
GS_VQE_res = GS_VQE.solve(problem)
```

```python id="-ziIfzxI6AIL"
# Resultats de la diagonalisation exacte
from qiskit_algorithms import NumPyMinimumEigensolver

numpy_solver = NumPyMinimumEigensolver()
NP_calc = GroundStateEigensolver(mapper, numpy_solver)
GS_NP_res = NP_calc.solve(problem)
```

<!-- #region id="fHfopi-S6AIL" -->
### Class `qiskit_nature_pyscf.PySCFGroundStateSolver`

Nous allons maintenant utiliser le plugin `qiskit_nature_pyscf` qui couple PySCF et Qiskit Nature.  C'est un solveur [FCI](https://en.wikipedia.org/wiki/Full_configuration_interaction) (Full Configuration Interaction) basé sur Qiskit-Nature qui permet à un utilisateur de PySCF (Python-based Simulations of Chemistry Framework) de tirer parti des algorithmes quantique implémentés dans Qiskit-Nature pour être utilisés à la place de leurs homologues classiques (dans un esprit similaire à l'intégration NWChemEx).

La classe `qiskit_nature_pyscf.PySCFGroundStateSolver` s'appuie sur le module ``fci`` de PySCF. Il n'utilise aucun algorithmes quantiques (puisqu'il les remplace dans le workflow de Qiskit-Naure) mais fournit à la place un utilitaire pour déboguer les workflows de calcul classique basés sur Qiskit-Nature.

Plus important encore, il fournit une implémentation plus efficace de ce que Qiskit-Nature réalise en utilisant la classe `qiskit_algorithms.NumPyMinimumEigensolver` en combinaison avec un ``filter_criterion``. Pour les états fondamentaux de spin autres que le singlet, l'utilisation des composants Qiskit-Nature est beaucoup plus complexe, alors que cette classe fournit une alternative facile à utiliser.
<!-- #endregion -->

```python id="C23NUR_A6AIL"
from pyscf import fci
try:
    from qiskit_nature_pyscf import PySCFGroundStateSolver
except:
    %pip install qiskit-nature-pyscf
    from qiskit_nature_pyscf import PySCFGroundStateSolver

au2kcalc = 627.509474
```

```python colab={"base_uri": "https://localhost:8080/"} id="1NNcJVfV6AIL" outputId="107c42ea-7fee-4e6f-cac0-5b49d7fb96d8"
fci_solver = fci.direct_uhf.FCI()
solver_pyscf = PySCFGroundStateSolver(fci_solver)

GS_FCI_res = solver_pyscf.solve(problem)
print(GS_FCI_res)
```

```python colab={"base_uri": "https://localhost:8080/"} id="rdIWCG776AIM" outputId="92d51f23-a072-4898-89d9-a3389f89e8e0"
print(f'error between the two exact calculations, {abs(GS_FCI_res.total_energies[0] - GS_NP_res.total_energies[0]) * au2kcalc} kcal/mol')
```

<!-- #region id="FD2Bj9Zv6AIM" -->
### Visualisation des résultats
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 112} id="hIMkCp236AIM" outputId="84c1bbe4-5c53-439a-d482-4dded4effa03"
error_VQE_NP = abs(GS_VQE_res.total_energies[0] - GS_NP_res.total_energies[0]) * au2kcalc
error_HF_NP = abs(GS_VQE_res.hartree_fock_energy - GS_NP_res.total_energies[0]) * au2kcalc

import pandas as pd

# Create the results dataframe
dict_res = {'FCI': [GS_FCI_res.total_energies[0], ''],
            'Numpy': [GS_NP_res.total_energies[0], ''],
            'VQE': [GS_VQE_res.total_energies[0], error_VQE_NP],
            'HF': [GS_VQE_res.hartree_fock_energy, error_HF_NP]}
df_GS = pd.DataFrame(dict_res, index = ['E_tot (a.u.)', 'Error (kcal/mol)'])

df_GS

```

<!-- #region id="39hNzUSF6AIM" -->
## Etats excités

Calculons maintenant les énergies des états excités de notre Hamiltonien moléculaire et en déduisons la
* la bande interdite ou l'écart entre le niveau fondamental $S_0$ et le premier niveau excité (le niveau triplet $T_1$),
* et et de l'énergie de fluorescence ou l'écart entre le niveau fondamental $S_0$ et le premier niveau excité singulet $S_1$.

![Molecule_HOMO-LUMO_diagram](Graphics/Molecule_HOMO-LUMO_diagram.png)

<!-- #endregion -->

<!-- #region id="Cqw_O1AH6AIN" -->
### Calcul avec `ExcitedStatesEigensolver`

Effectuons les calculs avec `NumPyEigensolver` avec le critère de filtre par défaut activé.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="E9iS1N8o6AIN" outputId="62c7fb53-d2a5-4c7b-dcf4-8737a87be1d0"
from qiskit_nature.second_q.algorithms import ExcitedStatesEigensolver
from qiskit_algorithms import NumPyEigensolver

numpy_solver = NumPyEigensolver(k=4, filter_criterion=problem.get_default_filter_criterion())

numpy_ES_solver = ExcitedStatesEigensolver(mapper, numpy_solver)
NP_ES = numpy_ES_solver.solve(problem)

print(NP_ES)
```

```python colab={"base_uri": "https://localhost:8080/"} id="mUw8M39_6AIN" outputId="6d342601-6d4f-4b79-cd21-270ff9d91bd6"
ET1_NP = NP_ES.total_energies[1]
ES1_NP = NP_ES.total_energies[2]
NP_gap = NP_ES.computed_energies[1] * au2ev
NP_f_energy = NP_ES.computed_energies[2] * au2ev

print(f'Total Numpy ES energy T1 = {ET1_NP} a.u.')
print(f'Total Numpy ES energy S1 = {ES1_NP} a.u.')
print(f'The bandgap obtained form Numpy ES calculations is : {NP_gap} eV')
print(f'The fluorescence energy obtained form Numpy ES calculations is : {NP_f_energy} eV')

```

```python colab={"base_uri": "https://localhost:8080/", "height": 175} id="3QiSS6SS6AIN" outputId="e2fd109e-5fa1-4b59-d80a-d728487790ba"
# Create the results dataframe
list_results_NP = [ET1_NP, ES1_NP, NP_f_energy, NP_gap]
dict_results_NP = {'Numpy ES': list_results_NP}
df_NP = pd.DataFrame(dict_results_NP,
                    index = ['ES energy T1 (a.u.)', 'ES energy S1 (a.u.)',
                            'f_energy ES1-ES0 (eV)', 'Gap ES1 - ET1 (eV)'])

df_NP
```

<!-- #region id="cn8oVwCZ6AIO" -->
### Calculs avec `PySCFGroundStateSolver`
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="l7FW49iC6AIO" outputId="88ebdc65-4a5e-4682-a5a8-4259e5af8f59"
fci_solver.nroots = 4
solver_fci = PySCFGroundStateSolver(fci_solver)
FCI_ES = solver_fci.solve(problem)
print(FCI_ES)
```

```python colab={"base_uri": "https://localhost:8080/"} id="kQqA_vNl6AIO" outputId="3a1d4989-2560-4395-c505-7eae97ce5851"
ET1_FCI = FCI_ES.total_energies[1]
ES1_FCI = FCI_ES.total_energies[2]
FCI_gap = FCI_ES.computed_energies[1] * au2ev
FCI_f_energy = FCI_ES.computed_energies[2] * au2ev

print(f'Total FCI ES energy T1 = {ET1_FCI} a.u.')
print(f'Total FCI ES energy S1 = {ES1_FCI} a.u.')
print(f'The bandgap obtained form FCI ES calculations is : {FCI_gap} eV')
print(f'The fluorescence energy obtained form FCI ES calculations is : {FCI_f_energy} eV')
```

```python colab={"base_uri": "https://localhost:8080/", "height": 175} id="psPxAsf86AIO" outputId="259ea253-7a4d-499b-9152-672bd4e4f138"
# Create the results dataframe
list_results_FCI = [ET1_FCI, ES1_FCI, FCI_f_energy, FCI_gap]
dict_results_FCI = {'FCI ES': list_results_FCI}
df_FCI = pd.DataFrame(dict_results_FCI,
                    index = ['ES energy T1 (a.u.)', 'ES energy S1 (a.u.)',
                            'f_energy ES1-ES0 (eV)', 'Gap ES1 - ET1 (eV)'])

df_FCI
```

<!-- #region id="kZoHQQQ-6AIP" -->
### Calculs avec `QEOM`

 Puisque nous avons déjà défini le système, nous avons besoin d'accéder à l'énergie d'excitation en utilisant l'[algorithme quantique d'équation du mouvement (qEOM)](https://arxiv.org/abs/1910.12890).

La classe `qiskit_nature.second_q.algorithms.QEOM` implémente cet algorithme qui approxime les propriétés de l'état excité d'un problème en utilisant des mesures supplémentaires sur l'état fondamental fournies par un objet `GroundStateSolver`. La précision de la méthode `GroundStateSolver.solve` pour l'approximation de l'état fondamental affecte directement la précision de l'algorithme qEOM pour le même problème. Les excitations sont utilisées pour construire un sous-espace linéaire dans lequel un problème de valeurs propres pour l'Hamiltonien projeté sera résolu. Cette méthode fonctionne généralement bien pour calculer les états excités les plus bas d'un problème. Les énergies des états excités sont calculées par défaut dans cet algorithme pour tous les états excités.
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/"} id="DXnIO0z96AIP" outputId="28571f6c-2f16-49e7-e01a-b6dff8f3c736"
from qiskit_nature.second_q.algorithms import QEOM
from qiskit_nature.second_q.algorithms.excited_states_solvers.qeom import EvaluationRule

qeom_ES_solver = QEOM(GS_VQE, Estimator(), "sd", EvaluationRule.ALL)
qeom_ES = qeom_ES_solver.solve(problem)
print(qeom_ES)
```

```python colab={"base_uri": "https://localhost:8080/"} id="5axFMtxs6AIP" outputId="00423e14-531d-40f8-fb67-6bea98ae11e3"
ET1 = qeom_ES.total_energies[1]
ES1 = qeom_ES.total_energies[2]
qeom_gap = qeom_ES.computed_energies[1] * au2ev
qeom_f_energy = qeom_ES.computed_energies[2] * au2ev

print(f'Total QEOM ES energy T1 = {ET1} a.u.')
print(f'Total QEOM ES energy S1 = {ES1} a.u.')
print(f'The bandgap obtained form Py ES calculations is : {qeom_gap} eV')
print(f'The fluorescence energy obtained form Py ES calculations is : {qeom_f_energy} eV')
```

```python colab={"base_uri": "https://localhost:8080/", "height": 175} id="GaDbauS66AIP" outputId="caf4ca1c-3e19-4434-d832-db9f45affbee"
# Create the results dataframe
list_results_QEOM = [ET1, ES1, qeom_f_energy, qeom_gap]
dict_results_QEOM = {'QEOM': list_results_QEOM}
df_QEOM = pd.DataFrame(dict_results_QEOM,
                    index = ['ES energy T1 (a.u.)', 'ES energy S1 (a.u.)',
                            'f_energy ES1-ES0 (eV)', 'Gap ES1 - ET1 (eV)'])

df_QEOM
```

<!-- #region id="skUd4Bnv6AIQ" -->
### Synthèse des résultats sur les états excités
<!-- #endregion -->

```python colab={"base_uri": "https://localhost:8080/", "height": 175} id="37ZoBuwM6AIQ" outputId="f3dc88e2-d3ba-4896-9cdf-f688cf0fed15"
# Visualisation of results

# Create the results dataframe
dict_res = {'FCI ES': list_results_FCI, 'Numpy ES': list_results_NP,
            'QEOM': list_results_QEOM
            }
df_ES = pd.DataFrame(dict_res,
                     index = ['ES energy T1 (a.u.)', 'ES energy S1 (a.u.)',
                            'f_energy ES1-ES0 (eV)', 'Gap ES1 - ET1 (eV)'])

df_ES

```

```python

```
