
project(san)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm -Wall -O2")
include_directories(ionic_currents)

set(SOURCE_FILES 
 ./ionic_currents/ca.c         
 ./ionic_currents/icat_t.c   
 ./ionic_currents/inaca.c	./ionic_currents/isus_t.c
 ./ionic_currents/ical_t.c	./ionic_currents/ifACH_t.c	./ionic_currents/ikr_t.c	./ionic_currents/ina_t.c  
 ./ionic_currents/clocks.c	./ionic_currents/icap.c		./ionic_currents/ip.c		./ionic_currents/ito_t.c
 ./ionic_currents/ibg.c		./ionic_currents/if_t.c		./ionic_currents/iks_t.c	./ionic_currents/pars.c
 ./ionic_currents/icalACH_t.c	./ionic_currents/ikACH_t.c  
 main3.c)
  
add_executable(san ${SOURCE_FILES})
