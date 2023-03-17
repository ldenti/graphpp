# graphpp - Graph PingPong
Specific strings computation from variation graphs.

### Prerequisites
* [vg](https://github.com/vgteam/vg)

### Compilation
```
git clone https://github.com/ldenti/graphpp.git
cd graphpp
mkdir build ; cd build
cmake ..
make
cd ..
```

### Usage
```
# index a graph with vg
vg index -g <graph.gcsa2> <graph.vg>
# compute specific strings
graphpp [-f FLANK] <graph.gcsa2> <query.fx>
```

### TODO
- [ ] example
- [ ] multithreading