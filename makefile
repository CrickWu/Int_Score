all:
	cd source_code && make
Release:
	cd source_code && make
	cp ./source_code/is_score.exe ./Int_score/bin/Release/Int_score.exe
clean:
	cd source_code && make clean
