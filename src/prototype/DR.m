globals;
[XX YY] = meshgrid(linspace(0,1,N+1),linspace(0,1,Q+1));
% Initialisation 
Zm = zeros(Q+1,N+1); Zf =zeros(Q+1,N+1);
Wm0 = zeros(Q+1,N+1); Wf0 =zeros(Q+1,N+1); Wm1 = zeros(Q+1,N+1); Wf1 =zeros(Q+1,N+1); Wm2 = zeros(Q+1,N+1); Wf2 = zeros(Q+1,N+1);
Zmbar = zeros(Q+1,N+2); Zfbar = zeros(Q+2,N+1); 
Wmbar0 = zeros(Q+1,N+2); Wfbar0 = zeros(Q+2,N+1); Wmbar1 = zeros(Q+1,N+2); Wfbar1 = zeros(Q+2,N+1); Wmbar2 = zeros(Q+1,N+2); Wfbar2 = zeros(Q+2,N+1); 
% Itérations
for i = 1:10
	[Zmbar,Zfbar,Zm,Zf] = proxG2(Wmbar0,Wfbar0,Wm0,Wf0);
	
	surf(XX,YY,Zf)
	xlabel('x')
	ylabel('t')
	drawnow
	[Wmbar1,Wfbar1,Wm1,Wf1] = proxG1(Wmbar0,Wfbar0,Wm0,Wf0,gamma);
	Wmbar1 = 2*Wmbar1 - Wmbar0; Wfbar1 = 2*Wfbar1 - Wfbar0;
	Wm1 = 2*Wm1 - Wm0; Wf1 = 2*Wf1 - Wf0;
	
	[Wmbar2,Wfbar2,Wm2,Wf2] = proxG2(Wmbar1,Wfbar1,Wm1,Wf1);
	Wmbar2 = 2*Wmbar2 - Wmbar1; Wfbar2 = 2*Wfbar2-Wfbar1; 
	Wm2 = 2*Wm2-Wm1; Wf2=2*Wf2-Wf1;
	
	Wmbar0 = (1-0.5*alpha)*Wmbar0 + 0.5*alpha*Wmbar2;
	Wfbar0 = (1-0.5*alpha)*Wfbar0 + 0.5*alpha*Wfbar2;
	Wm0 = (1-0.5*alpha)*Wm0 + 0.5*alpha*Wm2;
	Wf0 = (1-0.5*alpha)*Wf0 + 0.5*alpha*Wf2;
end