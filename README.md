# The Fungi on Living Plant Tissue, Fungi on Plant, and Locally-Co-occurring Fungi databases
<br>

## Description
### Three databases containing unique instances of co-occurrence among plant and fungal species, according to available biodiversity records. 
* The __Fungi on Living Plant Tissue__ database ([csv.zip](/Final_comb_DBs/FUNGI.ON.LIVING.PLANT.csv.zip) or [RDS](/Final_comb_DBs/FUNGI.ON.LIVING.PLANT.rds)) only includes occurrences of fungi _living directly on or within ___living___ plant tissue .<br><br>
* The __Fungi on Plant__ database ([csv.zip](/Final_comb_DBs/FUNGI.ON.PLANT.csv.zip) or [RDS](/Final_comb_DBs/FUNGI.ON.PLANT.rds)) is a more permissive database and includes direct contact for fungi on ___living or dead___ plant tissue. <br><br>
* The __Locally Co-Occurring Fungi__ database ([csv.zip](/Final_comb_DBs/LOCALLY.CO-OCCURRING.FUNGI.csv.zip) or [RDS](/Final_comb_DBs/LOCALLY.CO-OCCURRING.FUNGI.rds)) is the most permissive database and includes records of plant and fungal species indicating local geographical co-occurrence, ___whether or not direct tissue contact between organisms occurs___. <br>
<br>

Database files available in [Final_comb_DBs](/Final_comb_DBs).
<br><br>

## Read databases directly into R:
#### Fungi on Living Plant Tissue database:
```
living_plant<-readRDS(url("https://raw.githubusercontent.com/nicholasbard/fungi-plant/main/Final_comb_DBs/FUNGI.ON.LIVING.PLANT.rds"))
```

#### Fungi on Plant database:
```
on_plant<-readRDS(url("https://raw.githubusercontent.com/nicholasbard/fungi-plant/main/Final_comb_DBs/FUNGI.ON.PLANT.rds"))
```

#### Locally Co-Occurring Fungi database
```
local_co-occur<-readRDS(url("https://raw.githubusercontent.com/nicholasbard/fungi-plant/main/Final_comb_DBs/LOCALLY.CO-OCCURRING.FUNGI.rds"))
```
<br><br>
## Additional items (not essential for database use):

#### Protocols used for __Fungi on Living Plant Tissue__ database development are included in the [Rscript](/Rscript) directory. 
- [Database.Living.Tissue.Clean.2025.01.R](/Rscript/Database.Living.Tissue.Clean.2025.01.R) is the final protocol used for __Fungi on Living Plant__ database preparation and uses the [Cleaned_data](/Cleaned_data) created by using the remaining protocols. Additional data files from external sources used during database creation are located in [ext_data](/ext_data).

#### Fungi On Plant database development [protocol](/On.DB/Database.On_plant.R) and Locally Co-occurring Fungi database development [protocol](/Local.DB/Database.Local_co-occur.R). 
- Both use intermediate files (available upon request) created during the development of the __Fungi on Living Plant__ database (see above) and are found in [On.DB](/On.DB) and [Local.DB](/Local.DB), respectively.

#### R scripts for database analyses are provided in [Analyses](/Analyses).
- Used for data analysis and figures in upcoming paper, publication pending as of October 2025.





