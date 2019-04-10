
cmake_minimum_required(VERSION 3.5.0)

set(SDK_DIR ${CMAKE_SOURCE_DIR}/../../../SDK)

set(SRC_DIR ${SDK_DIR}/src)

set(TOOLCHAIN_PREFIX arm-none-eabi-)
set(CMAKE_TOOLCHAIN_FILE ${SRC_DIR}/arm-gcc-toolchain.cmake)

project(main.elf)

# Set a default build type if none was specified
set(default_build_type "Release")

# We need c99
set(CMAKE_C_FLAGS "-std=c99")

# options


# Compile flags
set(BFLAGS "-mcpu=cortex-m4 -mthumb -mlittle-endian -mfloat-abi=hard -mfpu=fpv4-sp-d16 -O1 -fmessage-length=0 -ffunction-sections -fdata-sections -ffreestanding -fsingle-precision-constant -Wall -Wextra -DARM_MATH_CM4 -D__FPU_PRESENT -Wl,--no-wchar-size-warning")
add_definitions(${BFLAGS} -DUSE_STDPERIPH_DRIVER -DHSE_VALUE=8000000 -DKODUINO_ISRS_INLINE -DUSING_RTOS -MMD -MP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fabi-version=0 -fno-exceptions -fno-rtti -fno-use-cxa-atexit -fno-threadsafe-statics")

set(CMAKE_EXE_LINKER_FLAGS "${BFLAGS} -fmessage-length=0 -fsigned-char -ffunction-sections -fdata-sections -ffreestanding -fno-move-loop-invariants -Wall -Wextra -nostartfiles -Xlinker --gc-sections --specs=nano.specs -u _printf_float")

# Include from Koduino directories
set(KODUINO_DIR ${SDK_DIR}/src/koduino/stm32)
include_directories(
	${SDK_DIR}/src
	${SDK_DIR}/inc
 	${KODUINO_DIR}/cores/arduino 
	${KODUINO_DIR}/system/CMSIS/Include 
 	${KODUINO_DIR}/system/trace
)

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
	add_definitions(-DUSE_FULL_ASSERT)
endif()

# Set F4
set(MCU_SERIES STM32F4xx)
set(MCU_PART STM32F446xx)
set(CLKSRC HSE_8M)
set(KODUINO_VARIANT f446)
set(LINKER_SCRIPT ${KODUINO_DIR}/variants/${KODUINO_VARIANT}/STM32F446ZETx_FLASH.ld)

set(STDPERIPH_SRC_DIR ${KODUINO_DIR}/system/${MCU_SERIES}_StdPeriph_Driver/src)
set(STDPERIPH_SRC_PREF ${STDPERIPH_SRC_DIR}/stm32f4xx_)
set(SRC_STDPERIPH
	${STDPERIPH_SRC_PREF}adc.c
	${STDPERIPH_SRC_PREF}dma.c
	${STDPERIPH_SRC_PREF}exti.c
	${STDPERIPH_SRC_PREF}gpio.c
	${STDPERIPH_SRC_PREF}misc.c
	${STDPERIPH_SRC_PREF}rcc.c
	${STDPERIPH_SRC_PREF}spi.c
	${STDPERIPH_SRC_PREF}syscfg.c
	${STDPERIPH_SRC_PREF}tim.c
	${STDPERIPH_SRC_PREF}usart.c
	${STDPERIPH_SRC_PREF}qspi.c
	${STDPERIPH_SRC_PREF}pwr.c
	)

# Create list of Koduino source files
file(GLOB SRC_KODUINO 
	 "${KODUINO_DIR}/cores/arduino/*.c" 
	 "${KODUINO_DIR}/cores/arduino/*.cpp"
	 "${KODUINO_DIR}/variants/${KODUINO_VARIANT}/*.c"
	 "${KODUINO_DIR}/variants/${KODUINO_VARIANT}/*.cpp"
	)
file(GLOB MCU_STARTUP 
	"${KODUINO_DIR}/variants/${KODUINO_VARIANT}/*.S"
	)

file(GLOB SRC_PROJECT 
"${CMAKE_SOURCE_DIR}/*.c" "${CMAKE_SOURCE_DIR}/*.cpp"
)

# Set MCU definitions
add_definitions(-DSERIES_${MCU_SERIES} -DCLKSRC_${CLKSRC} -D${MCU_PART} -DHSE_VALUE=8000000 -DNO_SDK)

# Include from Koduino paths
include_directories(
	${KODUINO_DIR}/system/${MCU_SERIES}_StdPeriph_Driver/inc 
	${KODUINO_DIR}/system/CMSIS/Device/ST/${MCU_SERIES}/Include 
	${KODUINO_DIR}/variants/${KODUINO_VARIANT}
	)

# Set assembly language on .S files
set_property(SOURCE ${MCU_STARTUP} PROPERTY LANGUAGE ASM)
enable_language(ASM)

# message(STATUS "SRC_ARCH " ${SRC_ARCH})

# Link math
link_directories(${KODUINO_DIR}/system/CMSIS/Lib/GCC)
link_libraries(arm_cortexM4lf_math)

# Set linker flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -T ${LINKER_SCRIPT}")

add_executable(main.elf ${SRC_PROJECT} ${SRC_KODUINO} ${SRC_STDPERIPH} ${MCU_STARTUP})

include(${SRC_DIR}/utils.cmake)
firmware_size(main.elf)

if (CMAKE_BUILD_TYPE MATCHES "Release")
	generate_object(main.elf .bin binary)
	
	# # Using ST-link
	# add_custom_command(TARGET main.elf POST_BUILD COMMAND "C:/Program Files (x86)/STMicroelectronics/STM32 ST-LINK Utility/ST-LINK Utility/ST-LINK_CLI.exe" -ME -P main.elf.bin 0x08000000 -Rst)

	# Using dfu
	# 1. Download dfu-util and add to path
	# 2. Use zadig to replace the driver for "STM32 BOOTLOADER" to WinUSB
	# add_custom_command(TARGET main.elf POST_BUILD COMMAND "touch" ../main.cpp)
	# add_custom_command(TARGET main.elf POST_BUILD COMMAND "${DFU_UTIL}" -a 0 -D main.elf.bin -s 0x08000000)
endif()

# 
