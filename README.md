## Our submission to the PACE 2024 Challenge to the Exact Track and the Parameterized Track
### Optil username: *mppeg*

Michael Jünger, Paul Jünger, Petra Mutzel, Gerhard Reinelt

June 9th, 2024

### Exact solver based on branch-and-cut

**Requirements:**
- Cbc Version 2.10.7
- CMake (we used version 3.22.1)

### To build the project, just run
```bash
cmake .
make
```

### Example usage
```bash
cat tiny_test_set/website_20.gr | ./oscm
# The solution will be printed to stdout
```

### Our results on the exact-public testset
The solver obtains provably optimal solutions to 99 of 100 instances in the required timelimit of 30 minutes.\
Tested on one thread of the _Intel Xeon "Sapphire Rapids" 48-core/96-threat 2.10GHz_ processor
|Instance | #Crossings|Time|
|---------|--------|---|
|exact-public/1.gr |1482| 0m 0s|
|exact-public/2.gr |3080| 0m 0s|
|exact-public/3.gr |6320| 0m 1s|
|exact-public/4.gr |6480| 0m 0s|
|exact-public/5.gr |9702| 0m 0s|
|exact-public/6.gr |11990| 0m 0s|
|exact-public/7.gr |14520| 0m 1s|
|exact-public/8.gr |16256| 0m 0s|
|exact-public/9.gr |19740| 0m 1s|
|exact-public/10.gr |16555| 0m 3s|
|exact-public/11.gr |20089| 0m 5s|
|exact-public/12.gr |829| 0m 0s|
|exact-public/13.gr |2744| 0m 0s|
|exact-public/14.gr |5316| 0m 0s|
|exact-public/15.gr |9500| 0m 0s|
|exact-public/16.gr |11068| 0m 1s|
|exact-public/17.gr |33251| 0m 12s|
|exact-public/18.gr |11841| 0m 0s|
|exact-public/19.gr |18104| 0m 0s|
|exact-public/20.gr |14897| 0m 0s|
|exact-public/21.gr |5176| 0m 0s|
|exact-public/22.gr |6777| 0m 0s|
|exact-public/23.gr |8590| 0m 0s|
|exact-public/24.gr |7686| 0m 0s|
|exact-public/25.gr |8139| 0m 0s|
|exact-public/26.gr |10879| 0m 0s|
|exact-public/27.gr |3230| 0m 0s|
|exact-public/28.gr |1559| 0m 0s|
|exact-public/29.gr |2776| 0m 0s|
|exact-public/30.gr |15024| 0m 0s|
|exact-public/31.gr |22312| 0m 0s|
|exact-public/32.gr |20873| 0m 0s|
|exact-public/33.gr |20724| 0m 0s|
|exact-public/34.gr |23408| 0m 0s|
|exact-public/35.gr |27740| 0m 0s|
|exact-public/36.gr |27022| 0m 0s|
|exact-public/37.gr |31948| 0m 0s|
|exact-public/38.gr |25208| 0m 1s|
|exact-public/39.gr |198926| 0m 1s|
|exact-public/40.gr |227764| 0m 0s|
|exact-public/41.gr |221630| 0m 0s|
|exact-public/42.gr |257869| 0m 1s|
|exact-public/43.gr |275442| 0m 1s|
|exact-public/44.gr |326396| 0m 1s|
|exact-public/45.gr |222924| 0m 0s|
|exact-public/46.gr |248405| 0m 0s|
|exact-public/47.gr |293935| 0m 1s|
|exact-public/48.gr |305888| 0m 1s|
|exact-public/49.gr |277023| 0m 0s|
|exact-public/50.gr |106802| 0m 3s|
|exact-public/51.gr |97850| 0m 1s|
|exact-public/52.gr |152556| 0m 1s|
|exact-public/53.gr |187314| 0m 3s|
|exact-public/54.gr |213217| 0m 5s|
|exact-public/55.gr |82205| 0m 0s|
|exact-public/56.gr |100013| 0m 0s|
|exact-public/57.gr |173013| 0m 0s|
|exact-public/58.gr |188442| 0m 9s|
|exact-public/59.gr |227475| 0m 3s|
|exact-public/60.gr |317024| 0m 3s|
|exact-public/61.gr |347582| 0m 3s|
|exact-public/62.gr |444898| 0m 9s|
|exact-public/63.gr |56563| 0m 8s|
|exact-public/64.gr |105838| 0m 7s|
|exact-public/65.gr |993019| 0m 1s|
|exact-public/66.gr |257876| 0m 3s|
|exact-public/67.gr |317718| 0m 5s|
|exact-public/68.gr |107438| 6m 22s|
|exact-public/69.gr |116996| 17m 17s|
|exact-public/70.gr |117037| 0m 0s|
|exact-public/71.gr |132493| 0m 0s|
|exact-public/72.gr |176033| 0m 0s|
|exact-public/73.gr |599603| 0m 26s|
|exact-public/74.gr |145468| 6m 20s|
|exact-public/75.gr |215824| 4m 20s|
|exact-public/76.gr |286207| 2m 37s|
|exact-public/77.gr |120099| 0m 18s|
|exact-public/78.gr |126862| 0m 19s|
|exact-public/79.gr |152071| 0m 36s|
|exact-public/80.gr |182715| 0m 34s|
|exact-public/81.gr |188778| 0m 43s|
|exact-public/82.gr |187569| 0m 32s|
|exact-public/83.gr |125099| 0m 0s|
|exact-public/84.gr |184166| 0m 0s|
|exact-public/85.gr |92759| 0m 0s|
|exact-public/86.gr |200617| 0m 0s|
|exact-public/87.gr |236782| 0m 0s|
|exact-public/88.gr |241803| 0m 0s|
|exact-public/89.gr |236418| 0m 0s|
|exact-public/90.gr |257813| 0m 0s|
|exact-public/91.gr |268908| 0m 0s|
|exact-public/92.gr |-| timeout|
|exact-public/93.gr |302803| 1m 57s|
|exact-public/94.gr |307447| 2m 37s|
|exact-public/95.gr |303429| 3m 30s|
|exact-public/96.gr |251921| 0m 51s|
|exact-public/97.gr |242361| 0m 0s|
|exact-public/98.gr |224831| 0m 1s|
|exact-public/99.gr |287587| 0m 0s|
|exact-public/100.gr |346841| 0m 0s|

_total (excluding timeout): 50m 55s_