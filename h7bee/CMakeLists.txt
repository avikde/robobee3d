
cmake_minimum_required(VERSION 3.5.0)

set(TOOLCHAIN_PREFIX arm-none-eabi-)
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/../arm-gcc-toolchain.cmake)

project(h7bee)

# Set a default build type if none was specified
set(default_build_type "Release")

# MCU options
set(MCU "-mcpu=cortex-m7 -mthumb -mfpu=fpv5-d16 -mfloat-abi=hard")
set(OPT -Og)

# defines
add_compile_definitions(USE_HAL_DRIVER STM32H743xx)
include_directories(
	${CMAKE_SOURCE_DIR}/Inc
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Inc
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Inc/Legacy
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/include
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/CMSIS_RTOS
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/portable/GCC/ARM_CM4F
	${CMAKE_SOURCE_DIR}/Drivers/CMSIS/Device/ST/STM32H7xx/Include
	${CMAKE_SOURCE_DIR}/Drivers/CMSIS/Include
)

# options
set(BFLAGS "${MCU} ${OPT} -Wall -fdata-sections -ffunction-sections -fsingle-precision-constant -MMD -MP")
set(CMAKE_C_FLAGS ${BFLAGS})
set(CMAKE_CXX_FLAGS ${BFLAGS})
if (CMAKE_BUILD_TYPE MATCHES "Debug")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -gdwarf-2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -gdwarf-2")
endif()

set(LDSCRIPT STM32H743VITx_FLASH.ld)
link_libraries(-lc -lm -lnosys)
set(CMAKE_EXE_LINKER_FLAGS "-specs=nano.specs  -u _printf_float -T${CMAKE_SOURCE_DIR}/${LDSCRIPT} -Wl,-Map=mb80.map,--cref -Wl,--gc-sections")

file(GLOB SRC_HAL
	${CMAKE_SOURCE_DIR}/Src/main.c
	${CMAKE_SOURCE_DIR}/Src/freertos.c
	${CMAKE_SOURCE_DIR}/Src/stm32h7xx_it.c
	${CMAKE_SOURCE_DIR}/Src/stm32h7xx_hal_msp.c
	${CMAKE_SOURCE_DIR}/Src/stm32h7xx_hal_timebase_tim.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_cortex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_spi.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_spi_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_tim.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_tim_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_uart.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_uart_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_rcc.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_rcc_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_flash.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_flash_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_gpio.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_hsem.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_dma.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_dma_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_mdma.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_pwr.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_pwr_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_i2c.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_i2c_ex.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_adc.c
	${CMAKE_SOURCE_DIR}/Drivers/STM32H7xx_HAL_Driver/Src/stm32h7xx_hal_adc_ex.c
	${CMAKE_SOURCE_DIR}/Src/system_stm32h7xx.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/croutine.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/event_groups.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/list.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/queue.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/stream_buffer.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/tasks.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/timers.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/CMSIS_RTOS/cmsis_os.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/portable/MemMang/heap_4.c
	${CMAKE_SOURCE_DIR}/Middlewares/Third_Party/FreeRTOS/Source/portable/GCC/ARM_CM4F/port.c
)
file(GLOB SRC_APP
	${CMAKE_SOURCE_DIR}/app/*.c
	${CMAKE_SOURCE_DIR}/app/*.cpp
)
set(MCU_STARTUP ${CMAKE_SOURCE_DIR}/startup_stm32h743xx.s)
# Set assembly language on .S files
set_property(SOURCE ${MCU_STARTUP} PROPERTY LANGUAGE ASM)
set_source_files_properties(${MCU_STARTUP} PROPERTIES COMPILE_FLAGS "-x assembler-with-cpp")
enable_language(ASM)

add_executable(h7bee.elf ${SRC_HAL} ${SRC_APP} ${SRC_SOEM} ${MCU_STARTUP})
# Post-build stuff
add_custom_command(TARGET h7bee.elf POST_BUILD 
COMMAND ${TOOLCHAIN_PREFIX}size h7bee.elf
COMMAND ${TOOLCHAIN_PREFIX}objcopy -O ihex h7bee.elf h7bee.hex
COMMAND ${TOOLCHAIN_PREFIX}objcopy -O binary -S h7bee.elf h7bee.bin
COMMENT "Post-build"
)

