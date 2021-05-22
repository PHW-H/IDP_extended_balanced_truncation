%Random model generator
function [Rl,Rc,Cc,Ll,ModT,M] = Random_model_generator(n,setMT,setRe,saveM)

Rl=randi([0 2000],n,1);
Rc=randi([0 2000],n,1);
Cc=randi([1 5000],n,1);
Cc=Cc*10^-6;
Ll=randi([50 15000],n,1);
Ll=Ll*10^-6;
if setMT==0
    ModT=randi([1 3]);
else
    ModT=setMT;
end
if setRe==0
    M=randi([0.1 0.5]);
else
    M=setRe;
end
if saveM==0
    save('Model_auto_save','Rl','Rc','Ll','Cc','ModT','M','n')
end
end
