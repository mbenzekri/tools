# OSMGRAF

Ce dépôt contient trois programmes principaux :

>**build_graph** construit un graphe routier orienté à partir d’un extrait OSM (PBF), en identifiant les nœuds d’intersection/terminaux, en découpant les ways en arcs, en calculant les longueurs, en gérant les sens uniques et en intégrant les restrictions de tournant, puis en écrivant un fichier binaire de graphe (format GRF2). Il peut aussi exporter des fichiers FlatGeobuf (nœuds, arcs, restrictions).

>**check_graph** lit ce fichier binaire et effectue des vérifications de cohérence + statistiques + connectivité.

>**search_graph** calcule tous les plus courts chemins entre un ensemble de noeuds d'interet (donnes via leurs node_id OSM), en utilisant Dijkstra + Radix Heap et en ecrivant les resultats en flux (texte ou binaire).

Les details ci-dessous s'appuient sur README.md, build_graph.cpp, check_graph.cpp et search_graph.cpp.

## Détails : build_graph

**But** : extraire un graphe routier à partir d’un fichier OSM PBF et produire un fichier binaire optimisé pour l’algorithme de routage.
Pipeline en plusieurs passes (résumé) :

1. Compter l’usage des nœuds des ways highway=* et repérer les extrémités.

2. Collecter les coordonnées et assigner un index aux nœuds de graphe (intersections/terminaux).

3. Découper les ways en arcs entre nœuds de graphe consécutifs, calculer la longueur (haversine), créer des arcs orientés en tenant compte des tags oneway et junction=roundabout.

4. Lire les restrictions de tournant (relations type=restriction, via = node), puis construire la liste des transitions interdites.

5. Construire le CSR (Compressed Sparse Row) et écrire le fichier binaire GRF2 qui contient : N, M, offsets, to, len, class, way_id, et des buffers pré-alloués pour dist/label/pred, puis les paires interdites.

**Optionnel** : export FlatGeobuf des nœuds/edges/restrictions via GDAL (--fgb).

**Sources**: build_graph.cpp (classes NodeUseCounter, NodeCollector, WayToArcs, RestrictionCollector, fonctions build_csr, build_forbidden_turns, write_graph_v2_with_restrictions, export_*, et main) et le pipeline décrit dans PROJECT.md. [build_graph.cpp] [PROJECT.md]

## Détails : check_graph

**But** : valider et analyser le graphe binaire produit.

Actions principales :

- Lit N et M, puis les tableaux CSR (offsets, to, len_cm, class, way_id).
- Vérifie la cohérence CSR : offsets[0]==0, offsets[N]==M, monotonie, to[e] < N.
- Calcule des statistiques : min/max/moyenne du degré sortant et des longueurs d’arcs.
- Calcule la connectivité en vue non orientée (BFS en utilisant arcs sortants + entrants) et donne nombre de composantes et taille de la plus grande.

**Source** : check_graph.cpp (fonction main, vérifications CSR, stats, BFS de connectivité). [check_graph.cpp]

## Details : search_graph

**But** : calculer les plus courts chemins entre D noeuds d'interet (node_id OSM) avec Dijkstra + Radix Heap, arret anticipe multi-cibles, et remise a zero rapide des tableaux.

Prerequis : une table de correspondance entre node_id OSM et node_id de graphe. Elle est lue :

Format du fichier nodes.txt (selection manuelle) :

- une ligne par point
- champs : node_id latitude longitude nom_du_lieu
- separateurs: espaces, ',' ou ';'
- le nom peut contenir des espaces (il est ignore par search_graph)
- soit via un fichier texte --node-map contenant des paires "graph_id osm_id" (ou "osm_id graph_id")
- soit via le fichier FlatGeobuf des noeuds (par defaut <base>_node.fgb, ou --node-fgb)

Usage (exemples) :

```
./search_graph out/france.graph in/nodes.txt out/paths.txt --text
./search_graph out/france.graph in/nodes.txt out/paths.bin --bin --node-fgb out/france_node.fgb
```

Format texte (une ligne par chemin) :

```
source_osm;destination_osm;distance_cm;node_id_1,node_id_2,...
```

Format binaire (par chemin) :

- int64 source_osm
- int64 destination_osm
- uint64 distance_cm (UINT64_MAX si pas de chemin)
- uint32 path_len
- int64[path_len] (liste des node_id OSM)

Limitations actuelles :

- pas de prise en compte des restrictions de tournant (les transitions interdites sont lues mais non appliquees au routage)

Entrées / sorties principales
**Entrée** : fichier OSM PBF (ex. fr-highway.osm.pbf).
**Sortie** : fichier binaire de graphe .graph en format GRF2.

**Optionnel** : fichiers FlatGeobuf (*_node.gfb, *_edge.gfb, *_restriction.gfb) si --fgb.

**Sources **: README.md, build_graph.cpp (commentaires et main). [README.md] [build_graph.cpp]

