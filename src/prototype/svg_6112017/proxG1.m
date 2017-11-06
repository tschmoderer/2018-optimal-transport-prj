function [Umbar,Ufbar,Vm,Vf] = proxG1(mbar,fbar,m,f,gamma)
	[Umbar Ufbar] = projC(mbar,fbar);
	[Vm Vf] = proxJ(m,f,gamma);
end