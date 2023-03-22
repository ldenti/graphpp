# graphpp - Graph PingPong
`graphpp` computes specific strings from variation graphs by implementing the [PingPong algorithm](https://academic.oup.com/bioinformaticsadvances/article/1/1/vbab005/6290003) on the [gcsa2 index](https://github.com/jltsiren/gcsa2).

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
./graphpp -h
./graphpp [-f FLANK] [-t THREADS] <graph.gcsa2> <query.fx>
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
- [X] ~avoid second assembly round if flank is 0~
- [X] do just one assembly round
- [X] support for fa and fq
- [X] multithreading
- [ ] improve multithreading
