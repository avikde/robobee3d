
cmake_minimum_required(VERSION 3.5.0)

set(TOOLCHAIN_PREFIX arm-none-eabi-)
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../arm-gcc-toolchain.cmake)

project(g4bee)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)

# UprightMPC2 controller
set(MPC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../template/uprightmpc2)
set(EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../SDK/thirdparty/eigen" CACHE STRING "Eigen location")
file(GLOB SRC_MPC ${MPC_DIR}/*.c)
# emosqp options
set(EMBEDDED 2)
set(PRINTING OFF CACHE BOOL "")
set(PROFILING OFF CACHE BOOL "")
set(DFLOAT ON CACHE BOOL "")
set(DLONG OFF CACHE BOOL "")
set(QDLDL_FLOAT_TYPE "float")
set(QDLDL_INT_TYPE "int")
set(QDLDL_BOOL_TYPE "int")
set(SRC_MPC_EIG ${MPC_DIR}/py/ematmult.cpp)
include_directories(${MPC_DIR} ${MPC_DIR}/py ${EIGEN_DIR})
add_library(uprightmpc2 STATIC ${SRC_MPC} ${SRC_MPC_EIG})

# Set a default build type if none was specified
set(default_build_type "Release")

# MCU options
set(MCU "-mcpu=cortex-m4 -mfpu=fpv4-sp-d16 -mfloat-abi=hard")
set(OPT -Og)

# defines
add_compile_definitions(USE_HAL_DRIVER STM32G431xx)
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/Inc
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Inc
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Inc/Legacy
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/CMSIS/Device/ST/STM32G4xx/Include
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/CMSIS/Include
)

# options
set(BFLAGS "${MCU} ${OPT} -Wall -fdata-sections -ffunction-sections -fsingle-precision-constant -MMD -MP")
set(CMAKE_C_FLAGS ${BFLAGS})
set(CMAKE_CXX_FLAGS ${BFLAGS})
if (CMAKE_BUILD_TYPE MATCHES "Debug")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -gdwarf-2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -gdwarf-2")
endif()

set(LDSCRIPT STM32G431KBTx_FLASH.ld)
link_libraries(-lc -lm -lnosys)
set(CMAKE_EXE_LINKER_FLAGS "-specs=nano.specs  -u _printf_float -T${CMAKE_CURRENT_SOURCE_DIR}/${LDSCRIPT} -Wl,--gc-sections")

set(SRC_HAL
	${CMAKE_CURRENT_SOURCE_DIR}/Src/main.c
	${CMAKE_CURRENT_SOURCE_DIR}/Src/stm32g4xx_it.c
	${CMAKE_CURRENT_SOURCE_DIR}/Src/stm32g4xx_hal_msp.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_gpio.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_tim.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_tim_ex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_pwr_ex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_uart.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_uart_ex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_rcc.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_rcc_ex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_flash.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_flash_ex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_flash_ramfunc.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_exti.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_dma.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_dma_ex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_pwr.c
	${CMAKE_CURRENT_SOURCE_DIR}/Drivers/STM32G4xx_HAL_Driver/Src/stm32g4xx_hal_cortex.c
	${CMAKE_CURRENT_SOURCE_DIR}/Src/system_stm32g4xx.c
)
set(SRC_APP
	${CMAKE_CURRENT_SOURCE_DIR}/app/_stdlib.c
	${CMAKE_CURRENT_SOURCE_DIR}/app/syscalls.c
	${CMAKE_CURRENT_SOURCE_DIR}/app/loop_update.cpp
)
set(MCU_STARTUP ${CMAKE_CURRENT_SOURCE_DIR}/startup_stm32g431xx.s)
# Set assembly language on .S files
set_property(SOURCE ${MCU_STARTUP} PROPERTY LANGUAGE ASM)
set_source_files_properties(${MCU_STARTUP} PROPERTIES COMPILE_FLAGS "-x assembler-with-cpp")
enable_language(ASM)

add_executable(g4bee.elf ${SRC_HAL} ${SRC_APP} ${MCU_STARTUP})
target_link_libraries(g4bee.elf uprightmpc2)
# Post-build stuff
add_custom_command(TARGET g4bee.elf POST_BUILD 
COMMAND ${TOOLCHAIN_PREFIX}size g4bee.elf
COMMAND ${TOOLCHAIN_PREFIX}objcopy -O ihex g4bee.elf g4bee.hex
COMMAND ${TOOLCHAIN_PREFIX}objcopy -O binary -S g4bee.elf g4bee.bin
COMMENT "Post-build"
)

# make -j16 && "C:\Program Files\STMicroelectronics\STM32Cube\STM32CubeProgrammer\bin\STM32_Programmer_CLI.exe" -c port=SWD -w g4bee.hex -hardRst