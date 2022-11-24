# 源文件路径为   ./src
# 目标文件路径为 ./obj
# 当前路径下创建目标文件路径文件夹
$(shell mkdir -p $(PACKAGE_PATH) ./obj)

FC     = gfortran    # 指定编译器
OMP    = -fopenmp    # 加openMP编译指令

# 源文件 及 目标文件
file   = ./src/commondata.f90 ./src/Dispersion.f90 ./src/pretreatment.f90 ./src/GMF.f90 ./src/PSO.f90 ./src/MainDrive.f90
obj    = ./obj/commondata.o   ./obj/Dispersion.o   ./obj/pretreatment.o   ./obj/GMF.o   ./obj/PSO.o   ./obj/MainDrive.o
# 指定可执行文件名称 
target = GMFinv.run      

# 链接指定目标文件路径内的所有.o文件，并生成可执行程序
$(target):$(obj)
	$(FC) $(obj) -static-libgfortran -o $@ $(OMP)

# 编译指定路径（./src）内的所有源文件，并将目标文件生成到指定路径（.obj）
./obj/%.o:./src/%.f90
	$(FC) $(OMP) -c $< -o $@


.PHONY:
clean:
	rm -rf ./obj/*.o ./obj/*.mod
	rm -rf *.o *.mod
	rm -rf $(target)