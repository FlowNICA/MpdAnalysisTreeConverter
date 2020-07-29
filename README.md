# AnalysisTree converter for MPD (NICA)

## Installation

- Export necessary environment variables
```bash
source <path>/mpdroot/build/config.sh
export ANALYSISTREE_INC=<path>/AnalysisTree/install-cxx11/include/
export ANALYSISTREE_LIB=<path>/AnalysisTree/install-cxx11/lib/
```

- Get the source code:
```bash
git clone https://devel.mephi.ru/PEParfenov/MpdAnalysisTreeConverter.git
cd MpdAnalysisTreeConverter/
mkdir build
cd build/
cmake ..
make
```

## Usage

```bash
./MpdDst2AnalysisTree -i input_mpddst.root -o output_AnalysisTree.root
```
