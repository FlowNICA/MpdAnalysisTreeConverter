#!/bin/bash

#---Option: if set, script stops (with logout) when encounters errors
#set -e

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



#---Install gcc ver. 8.5.0
if [ ! -d "${INSTALLATION_PATH}/gcc-8.5.0" ]
then
  echo "Installing gcc ver. 8.5.0 ..."
  cd $INSTALLATION_PATH
  wget http://mirror.linux-ia64.org/gnu/gcc/releases/gcc-8.5.0/gcc-8.5.0.tar.gz
  tar xvfz gcc-8.5.0.tar.gz
  cd gcc-8.5.0/
  ./contrib/download_prerequisites
  ./configure --disable-multilib --prefix=${INSTALLATION_PATH}/gcc-8.5.0
  make -j${NPROC}
  make install
else
  echo "Directory gcc-8.5.0 already exists in the installation path. Skipping this step."
fi



#---Install boost ver. 1.77.0
if [ ! -d "${INSTALLATION_PATH}/boost_1_77_0" ]
then
  echo "Installing boost ver. 1.77.0 ..."
  cd $INSTALLATION_PATH
  wget https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.gz
  tar xvfz boost_1_77_0.tar.gz
  cd boost_1_77_0/
	export PATH=${INSTALLATION_PATH}/gcc-8.5.0/bin:$PATH
	export LD_LIBRARY_PATH=${INSTALLATION_PATH}/gcc-8.5.0/lib64/:$LD_LIBRARY_PATH
	export BOOST_ROOT=${INSTALLATION_PATH}/boost_1_77_0/
	./bootstrap.sh --prefix=${INSTALLATION_PATH}/boost_1_77_0/
	./b2 install
else
  echo "Directory boost_1_77_0 already exists in the installation path. Skipping this step."
fi



#---Install cmake ver. 3.18
if [ ! -d "${INSTALLATION_PATH}/cmake-3.18.0" ]
then
  echo "Installing cmake ver. 3.18.0 ..."
  cd $INSTALLATION_PATH
  wget https://github.com/Kitware/CMake/releases/download/v3.18.0/cmake-3.18.0.tar.gz
  tar xfvz cmake-3.18.0.tar.gz
  cd cmake-3.18.0/
  ./bootstrap --prefix=${INSTALLATION_PATH}/cmake-3.18.0
  make -j${NPROC}
  make install
  export CMAKE_BIN=${INSTALLATION_PATH}/cmake-3.18.0/bin/cmake
else
  echo "Directory cmake-3.18.0 already exists in the installation path. Skipping this step."
  export CMAKE_BIN=${INSTALLATION_PATH}/cmake-3.18.0/bin/cmake
fi



#---Install root 6.20 with cxx-17
if [ ! -d "${INSTALLATION_PATH}/root/build-cxx17" ]
then
  echo "Installing root ver. 6.20 with c++17 standard ..."
  cd $INSTALLATION_PATH
  mkdir root/
  cd root/
  git clone https://github.com/root-project/root.git src/
  cd src/
  git checkout v6-20-00-patches
  cd ../
  mkdir build-cxx17 install-cxx17
  cd build-cxx17/
  export CC=${INSTALLATION_PATH}/gcc-8.5.0/bin/gcc
  export CXX=${INSTALLATION_PATH}/gcc-8.5.0/bin/g++
  export LD_LIBRARY_PATH=${INSTALLATION_PATH}/gcc-8.5.0/lib64:$LD_LIBRARY_PATH
  $CMAKE_BIN -Dbuiltin_gsl=ON -Dmathmore=ON -Dbuiltin_xrootd=ON -Dbuiltin_glew=ON -Dbuiltin_afterimage=ON -Dbuiltin_vdt=ON -DCMAKE_INSTALL_PREFIX=${INSTALLATION_PATH}/root/install-cxx17/ -DCMAKE_CXX_STANDARD=17 ../src/
  make -j${NPROC} install
else
  echo "Directory root/build-cxx17 already exists in the installation path. Skipping this step."
fi



#---Install AnalysisTree core part (with c++17 standard)
if [ ! -d "${INSTALLATION_PATH}/AnalysisTree/build-cxx17" ]
then
  echo "Installing AnalysisTree with c++17 standard ..."
  cd $INSTALLATION_PATH
  if [ ! -d "${INSTALLATION_PATH}/AnalysisTree/" ]
  then
    git clone https://github.com/HeavyIonAnalysis/AnalysisTree.git
  fi
  cd AnalysisTree/
  git checkout release-1.0
  mkdir build-cxx17 install-cxx17
  cd build-cxx17/
  source ${INSTALLATION_PATH}/root/build-cxx17/bin/thisroot.sh
  export CC=${INSTALLATION_PATH}/gcc-8.5.0/bin/gcc
  export CXX=${INSTALLATION_PATH}/gcc-8.5.0/bin/g++
  export LD_LIBRARY_PATH=${INSTALLATION_PATH}/gcc-8.5.0/lib64:$LD_LIBRARY_PATH
  export CPATH=${INSTALLATION_PATH}/root/build/include/
  $CMAKE_BIN -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=${INSTALLATION_PATH}/AnalysisTree/install-cxx17/ ../
  make -j${NPROC} install
else
  echo "Directory AnalysisTree/build-cxx17 already exists in the installation path. Skipping this step."
fi



#---Install QnAnalysis with cxx17, root-6.20, gcc-8.5.0, boost-1.77.0
if [ ! -d "${INSTALLATION_PATH}/QnAnalysis" ]
then
  echo "Installing QnAnalysis (with cxx17 standard, root-6.20, gcc-8.5.0, boost-1.77.0) ..."
  cd $INSTALLATION_PATH
  git clone git@github.com:HeavyIonAnalysis/QnAnalysis.git
  cd QnAnalysis/
  mkdir build install
  cd build/
  source ${INSTALLATION_PATH}/root/build-cxx17/bin/thisroot.sh
  export PATH=${INSTALLATION_PATH}/gcc-8.5.0/bin:$PATH
	export LD_LIBRARY_PATH=${INSTALLATION_PATH}/gcc-8.5.0/lib64/:$LD_LIBRARY_PATH
	export BOOST_ROOT=${INSTALLATION_PATH}/boost_1_77_0/
  $CMAKE_BIN -DCMAKE_INSTALL_PREFIX=${INSTALLATION_PATH}/QnAnalysis/install -Dyaml-cpp_BUNDLED=ON ../
	make -j${NPROC} install
else
  echo "Directory QnAnalysis already exists in the installation path. Skipping this step."
fi



cd $STARTING_PWD
