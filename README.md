# FBMC met gediversifieerde minRAM

Omschrijving van de bestanden:
* **main.jl** bevat een market coupling model **m** en congestion management model **c**
* **data-basic** en **data-advanced** bevatten alle nodige datasets. Wisselen van datasets gaat gemakkelijk in main.jl
* de originele datasets zijn **df~, matrix~** en **data~**
* **sortData.py** bewerkt de datasets
* de bewerkte datasets zijn **branch, bus, incidence, load, plant, susceptance**
* andere **csv** files zijn de output van het market coupling en congestion management model
* **plotResults.py** geeft de resultaten weer in grafieken
