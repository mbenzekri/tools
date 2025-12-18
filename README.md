# Decoupage du graphe routier OSM en niveau d'importance

L'idée est d'affecter un niveau d'importance à chaque troncon routier OSM.
Les niveaux ainsi constitués doivent permettre

- des affichages selon l'échelle par exemple
- du calcul d'itinéraire sur différent niveau

Niveau:

0. national
1. régional
2. départemental
3. communal
4. total

- Les niveaux doivent consituer un graphe quasi-connexe
- Chaque niveau est imbriqué dans le suivant

WARNING:
Le graphe routier n'existe pas en tant que tel dans OSM.
Les routes sont des objets linéaire (way/way_id) composés de noeuds/point (node/node_id).
Les way (suite de node) peuvent contenir des intersections routière sur des nodes intermédiaire (hors extremité)

À ce stade :

- **pas de restrictions de tournant**
- **pas encore de gestion des sens uniques**
- priorité donnée à la **topologie correcte** et à la **qualité structurelle du graphe**

## Processus

- extraire le graphe routier depuis une archive OSM
    ``osmium tags-filter in/france-latest.osm.pbf w/highway -o out/fr-highway.osm.pbf``
- construire le graphe routier depuis l'extrait
    ``./tools/./build_graph out/fr-highway.osm.pbf out/france.graph``
- valider le graphe routier depuis l'extrait
    ``./tools/check_graph out/france.graph``
- générer les niveaux
    ``./tools/classify_graph ...``

## Données d’entrée

- Source : **OpenStreetMap**
- Fournisseur : **Geofabrik**
- Exemple utilisé pour les tests : **Alsace** (`alsace-latest.osm.pbf`) puis filtre `highway=*`

---

## Pipeline mis en place

### 1) Lecture OSM (.osm.pbf)

- Langage : **C++17**
- Bibliothèque : **libosmium**
- Lecture directe du format **PBF**

### 2) Définition des nœuds du graphe

Un **nœud de graphe** est défini comme :

- un nœud OSM **partagé par au moins deux ways routiers** (intersection)
- **ou** une **extrémité** de way routier

Tous les autres nœuds OSM (géométrie intermédiaire) sont utilisés **uniquement** pour le calcul des longueurs, mais **ne deviennent pas des nœuds du graphe**.

### 3) Découpage des ways en arcs

Pour chaque `way` portant le tag `highway=*` :

- parcours séquentiel de la liste des nodes OSM
- **découpage** entre deux nœuds de graphe consécutifs
- calcul de la longueur du tronçon = somme des distances (haversine) entre nodes OSM successifs
- génération de **2 arcs dirigés** (double sens, provisoire)

À ce stade :

- pas de prise en compte de `oneway`
- pas de filtrage `access`, `motor_vehicle`, etc.

### 4) Structure du graphe (CSR)

Le graphe est stocké en **CSR (Compressed Sparse Row)** orientée :

- `offsets[N+1]` : début des arcs sortants de chaque nœud
- `to[M]` : destination de chaque arc
- `len_cm[M]` : longueur en centimètres
- `class[M]` : classe de route (actuellement 0)
- `way_id[M]` : identifiant OSM du way d’origine

### 5) Fichier binaire produit

Le programme `build_graph.cpp` génère un fichier binaire auto-contenu :

```txt
[uint32] N
[uint32] M

[offsets   : uint32 × (N+1)]
[to        : uint32 × M]
[len_cm    : uint32 × M]
[class     : uint8  × M]
[way_id    : int64  × M]

[dist      : uint64 × N]  // initialisé à INF
[label     : uint32 × N]  // initialisé à NO_NODE
[pred      : uint32 × N]  // initialisé à NO_NODE
```

Les tableaux `dist / label / pred` sont **pré-alloués et initialisés**, afin que la phase Moore–Dijkstra puisse s’exécuter **sans réallocation**.

### 6) Audit du graphe

Un second programme `check_graph.cpp` permet de **valider** le graphe produit.

Vérifications effectuées :

- `offsets[0] == 0`
- `offsets[N] == M`
- monotonie de `offsets`
- `to[e] < N` pour tous les arcs
- degré sortant min / max / moyen
- longueurs min / max / moyenne
- nombre de segments de longueur nulle
- calcul des **composantes connexes** sur la vue non orientée (sortants + entrants)
- taille de la plus grande composante (GCC)
- nombre total de composantes

---

## Ce qui est volontairement absent (pour l’instant)

- restrictions de tournant (`type=restriction`)
- gestion des sens uniques (`oneway`, `junction=roundabout`)
- filtrage fin par accessibilité (`motor_vehicle`, `service`, etc.)
- routage effectif (Moore–Dijkstra non encore implémenté)

---

## État actuel

Le projet dispose maintenant :

- d’un **build de graphe fonctionnel** et reproductible
- d’un **format binaire stable**
- d’un **outil d’audit** (structure + stats + connectivité)
- de métriques objectives sur la qualité du graphe

La base est prête pour :

- extraction de la plus grande composante
- ajout progressif des règles routières (`oneway`, access…)
- implémentation du Moore–Dijkstra multi-source
