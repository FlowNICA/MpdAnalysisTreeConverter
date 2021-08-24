#!/bin/bash

#---Option: if set, script stops (with logout) when encounters errors
#set -e

export MPDROOT_CONFIG=/scratch2/parfenov/Soft/MPDROOT/mpdroot_kfparticles/build/config.sh

NPROC=`nproc`
echo "Number of available CPUs: $NPROC"

export STARTING_PWD=$PWD
export INSTALLATION_PATH=$PWD/testQn
mkdir -p $INSTALLATION_PATH

#---Save original environmant
export PATH_Qn_original=$PATH
export LD_LIBRARY_PATH_Qn_original=$LD_LIBRARY_PATH
export BOOST_ROOT_Qn_original=$BOOST_ROOT
export CC_Qn_original=$CC
export CXX_Qn_original=$CXX




#---Install root 6.20 with cxx-11
if [ ! -d "${INSTALLATION_PATH}/root/build-cxx11" ]
then
  echo "Installing root ver. 6.20 with c++11 standard ..."
  cd $INSTALLATION_PATH
  mkdir root/
  cd root/
  git clone https://github.com/root-project/root.git src/
  cd src/
  git checkout v6-20-00-patches
  cd ../
  mkdir build-cxx11 install-cxx11
  cd build-cxx11/
  source $MPDROOT_CONFIG
  # export CC=${INSTALLATION_PATH}/gcc-8.5.0/bin/gcc
  # export CXX=${INSTALLATION_PATH}/gcc-8.5.0/bin/g++
  # export LD_LIBRARY_PATH=${INSTALLATION_PATH}/gcc-8.5.0/lib64:$LD_LIBRARY_PATH
  cmake -Dbuiltin_gsl=ON -Dmathmore=ON -Dbuiltin_xrootd=ON -Dbuiltin_glew=ON -Dbuiltin_afterimage=ON -Dbuiltin_vdt=ON -DCMAKE_INSTALL_PREFIX=${INSTALLATION_PATH}/root/install-cxx11/ -DCMAKE_CXX_STANDARD=11 ../src/
  make -j${NPROC} install
else
  echo "Directory root/build-cxx11 already exists in the installation path. Skipping this step."
fi



#---Install AnalysisTree core part (with c++11 standard)
if [ ! -d "${INSTALLATION_PATH}/AnalysisTree/build-cxx11" ]
then
  echo "Installing AnalysisTree with c++11 standard ..."
  cd $INSTALLATION_PATH
  git clone https://github.com/HeavyIonAnalysis/AnalysisTree.git
  cd AnalysisTree/
  git checkout release-1.0
  mkdir build-cxx11 install-cxx11
  cd build-cxx11/
  export CC=$CC_Qn_original
  export CXX=$CXX_Qn_original
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_Qn_original
  source $MPDROOT_CONFIG
  # source ${INSTALLATION_PATH}/root/build-cxx11/bin/thisroot.sh
  cmake -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=${INSTALLATION_PATH}/AnalysisTree/install-cxx11/ ../
  make -j${NPROC} install
  export AnalysisTreeBase_LIBRARY_DIR=${INSTALLATION_PATH}/AnalysisTree/install-cxx11/lib
  export AnalysisTree_INCLUDE_DIR=${INSTALLATION_PATH}/AnalysisTree/install-cxx11/include
else
  echo "Directory AnalysisTree/build-cxx11 already exists in the installation path. Skipping this step."
  export AnalysisTreeBase_LIBRARY_DIR=${INSTALLATION_PATH}/AnalysisTree/install-cxx11/lib
  export AnalysisTree_INCLUDE_DIR=${INSTALLATION_PATH}/AnalysisTree/install-cxx11/include
fi



cd $STARTING_PWD
