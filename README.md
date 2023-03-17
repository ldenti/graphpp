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
./graphpp [-f FLANK] <graph.gcsa2> <query.fx>
```

### Example
```
cd example
vg convert -g graph.gfa -v > graph.vg
vg index -g graph.gcsa2 graph.vg
../graphpp graph.gcsa2 seqs.fa > specific.fa
```

### TODO
- [X] example
- [ ] avoid second assembly round if flank is 0
- [ ] multithreading