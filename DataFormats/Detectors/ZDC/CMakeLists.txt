# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(DataFormatsZDC
               SOURCES src/ChannelData.cxx src/BCData.cxx src/RecEvent.cxx
	       	       src/OrbitRawData.cxx src/OrbitRecData.cxx 
               PUBLIC_LINK_LIBRARIES O2::CommonConstants O2::CommonDataFormat
	       			     O2::ZDCBase ROOT::MathCore FairRoot::Base
				     O2::MathUtils ms_gsl::ms_gsl)

o2_target_root_dictionary(DataFormatsZDC
               HEADERS include/DataFormatsZDC/BCData.h include/DataFormatsZDC/ChannelData.h
	       	       include/DataFormatsZDC/RecEvent.h include/DataFormatsZDC/OrbitRawData.h
		       include/DataFormatsZDC/OrbitRecData.h)
