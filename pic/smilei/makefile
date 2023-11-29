default: setup test sim

setup:
	@echo "Setup"
	@mkdir "run/$(VAR)" -p

.PHONY: test
test:
	@echo "Test"
	@cd "run/$(VAR)" && mkdir "log" -p && \
	$(spack_load) && \
	$(HOME)/Smilei/smilei_test test.py > log/test.log

sim:
	@echo "Simlutation"
	@cd $(VAR) && mkdir "log" -p && \
	$(spack_load) && \
	$(HOME)/Smilei/smilei test.py > log/sim.log

clean:
	@echo "Cleanning"
	cd $(VAR) && rm -rf Fields*.h5  ParticleBinning*.h5  patch_load.txt  Probes*.h5  profil.txt  scalars.txt  smilei.py

job:
	qsub zijin.pbs

spack_load = spack env activate smilei-oneapi