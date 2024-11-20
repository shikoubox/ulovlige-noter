math2 := module()
description "Some useful tools for Math 2 exams";
option package;
export prik, kryds, len, vop, phase_portrait, conv_check, decomp, conv_radi, leibniz_criterion, tranf_func, sys_tranf, equiv, outarg, method2;
    	with(LinearAlgebra);
    	with(plots);
	with(plottools);
	with(ArrayTools);
    	with(DEtools);
    	
# 	interface(imaginaryunit=i);


    	prik:= proc(x::Vector,y::Vector)
    	description "prikprodukt mellem vektorer";
        	return VectorCalculus[DotProduct](x,y);
    	end proc;

    	kryds:= proc(x::Vector, y::Vector)
    	description "krydsprodukt mellem vektorer";
        	return convert(VectorCalculus[CrossProduct](x,y),Vector);
    	end proc;

    	len:= proc(a::{Vector})
	description "sqrt(v²)";
        	return sqrt(prik(a,a));
    	end proc;

    	vop:=proc(x)
    	description "vop";
        	op(convert(x,list))
    	end proc:

    	#Konvergens tjek af en række
    	#kunne være sejt at indføre ækvivalens kriteriet LOL
    	conv_check := proc(a)
    		local kvo;
    		#n'te-ledskriteriet
    		if limit(a, n = infinity) <> 0 then
    			print("Nej, fra n'te ledskriteriet sætn. 4.7, rækken er divergent");
		return ;
		else
			print("n'te-ledskriteriet siger ikke noget");
		end if;
		
		#Kvotientkriteriet
		kvo := subs(n = n + 1, abs(a))/abs(a);
		if limit(kvo, n = infinity) < 1 then
			print("Absolut konvergent ifølge kvotientkriteriet sætn. 4.30");
		return ;
		else
			print("Kvotientkriteriet siger ikke noget");
		end if;

		#Integralkriteriet
		if abs(evalf(int(abs(a), n = 1 .. infinity))) < infinity then
			print("Absolut konvergent iflg. integralkriteriet sætn. 4.33 og def. 4.26");
		elif abs(evalf(int(abs(a), n = 2 .. infinity))) = Float(infinity) then
			print("Betinget konvergent iflg. integralkriteriet sætn. 4.33");
		return ;
		else
			print("Divergent iflg. Integralkriteriet sætn. 4.33");
		return ;
		end if;
	
	end proc;

	#Leibniz criterion, is used for alternating series.
	leibniz_criterion := proc(b)
		if limit(abs(b),n=infinity) = 0 then
			if min(b)>0 then
				if b < subs(n=n+1,b) then
					print("konvergent");
				else print("Aftager ikke monotont");
				end if;
			else
				print("Har ikke kun positive led");
			end if;	
   		print("Konvergent ifølge leibniz' kriterie");
	 	else
		 	printf("Divergent, da det er en alternerende række, hvor det n'te led ikke konvergerer mod 0 men mod %a", limit(abs(b),n=infinity));
		end if;
	end proc:

	#Dekomposition
	decomp := proc(b)
		return convert(b,parfrac,x);
	end proc:

	#Faseportræt omkring (0,0)
	phase_portrait := proc(eqn1,eqn2)
		local ic1, ic2;
		ic1 := x1(0)=-0.5, x2(0)=1;ic2 := x1(0)=0.5, x2(0)=-0.5;
		DEplot([eqn1,eqn2],[x1(t),x2(t)], t=-1..1,[[ic1],[ic2]],stepsize=.05);
	end proc:
	
	#Metode (i)
	

	#Metode (ii) 
	#Approksimation jeg ikke ved hvornår man skal bruge (approksimation af fejl)
	method2 := proc(a,epsilon)
		print(ceil(max(solve(abs(subs(n=N+1,a)) = epsilon,N))));
	end proc:

	#Ækvivalent række tjek equiv(a, b) hvor a og b er de to rækker der skal sammenlignes.
	equiv := proc(a, b)
		if is(abs(limit(a/b, n = infinity)) < infinity) = true and 0 < abs(limit(a/b, n = infinity))
		then print("ækvivalente rækker");
		else print("ikke ækvivalente rækker");
		end if;
	end proc:

	#Konvergens radius fra kvotientkriteriet jvf. 4.30
	conv_radi := proc(a)
		local kvotient, li, rho;
		kvotient := abs(subs(n = n + 1, a)/a);
		print(kvotient);
		li := limit(kvotient, n = infinity);
		print(li);
		rho := solve(li = 1, x);
		print("radius of convergence is jvf. sætning 4.30. then "*'rho' = rho);
	end proc;

	#outarg ???? ved ikke hvad er
	outarg := proc(a)
		local arg1, arg2, arg3, argsin, argcos, argexp;
		arg1 := convert(indets(a, 'cos(anything)'), list);
		argcos := map2(op, 1, arg1);
		arg2 := convert(indets(a, 'sin(anything)'), list);
		argsin := map2(op, 1, arg2);
		arg3 := convert(indets(a, 'exp(anything)'), list);
		argexp := map2(op, 1, arg3);
		if numelems(arg1) = 0 and numelems(arg3) = 0 then
			argsin := argsin[];
			argcos := 0;
			argexp := 0;
		elif numelems(arg2) = 0 and numelems(arg3) = 0 then
			argcos := argcos[];
			argsin := 0;
			argexp := 0;
		end if;
		if numelems(arg1) = 0 and numelems(arg2) = 0 then 
			argexp := argexp[];
			argsin := 0;
			argcos := 0;
		elif type(argexp, list) = true and numelems(arg1) = 0 then
			argexp := argexp[];
			argsin := argsin[];
			argcos := 0;
		elif type(argexp, list) = true and numelems(arg2) = 0 then
			argexp := argexp[];
			argcos := argcos[];
			argsin := 0;
		end if;
		return [argsin, argcos, argexp];
	end proc;

	# tranf_func(a,b,d) hvor
	# 'a' er en række vektor indeholdene koefficienterne foran venstreside af ligningen (y'erne), 
	# "b" er en række vektor indeholdene koefficienter på højresiden (u'erne)
	# *d* er den påvirkning som systemet påtrykkes, Ønsker du kun at se overføringsfunktionen sættes d til 0!

	tranf_func := proc(a, b, d)
		local y, vec, vec1, vec2, u, H, lambda, argu, argsin, argcos, argexp, f, losolve, g, y1, c;
		y := exp(s*t);
		for f to numelems(a) - 1 do
			y := <diff(exp(s*t), t $ f), y>;
		end do;
		vec := a . y;
		if numelems(a) - 1 = 1 then
			vec1 := vec;
		else vec1 := vec[1];
		end if;
		
		u := exp(s*t);
		for f to numelems(b) - 1 do
			u := <diff(exp(s*t), t $ f), u>;
		end do;
		vec := b . u;
		if numelems(b) - 1 = 1 then
			vec2 := vec;
		elif numelems(b) - 1 = 0 then
			vec2 := u;
		else
			vec2 := vec[1];
		end if;
		
		H := s -> simplify(vec2/vec1);
		print("Overføringsfunktionen bliver da");
		print('H(s)' = H(s));
		print("s skal være forskellig fra de nedenstående værdier");
		losolve := solve(denom(H(s)));
		for f to degree(denom(H(s))) do
			print(losolve <> simplify([f], complex));
		end do;
		if d = 0 then return ;
		end if;
		
		argu := outarg(d);
		argsin := argu[1];
		argcos := argu[2];
		argexp := argu[3];
		
		print("den partikulære løsning bliver da");
		if d = sin(argsin) then
			lambda := argsin*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Im(subs(s = lambda/t, H(s))*exp(lambda)));
		elif d = cos(argcos) then
			lambda := argcos*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Re(subs(s = lambda/t, H(s))*exp(lambda)));
		elif d = exp(argexp) then
			lambda := argexp;
#			print("den partikulære løsning bliver da");
			y := subs(s = lambda/t, H(s))*exp(lambda);
		elif d = cos(argcos)*exp(argexp) then
			lambda := argexp + argcos*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Re(subs(s = lambda/t, H(s))*exp(lambda)));
		elif d = sin(argsin)*exp(argexp) then
			lambda := argexp + argsin*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Im(subs(s = lambda/t, H(s))*exp(lambda))); 
		end if;
	end proc;

	# For et system af diff. lign. bruges funktionen sys_tranf.
	sys_tranf := proc(a, b, c, d)
		local num, funk, H, argcos, argsin, argexp, y, lambda, argu;
		num := Size(a, 1);
		funk := simplify(-(((Transpose(c)) . (1/(a - s*IdentityMatrix(num)))) . b));
		H := s -> funk; print("overførings funktionen bliver da"*'H(s)' = H(s));
		if d = 0 then return ; end if;
		argu := outarg(d);
		argsin := argu[1];
		argcos := argu[2]; 
		argexp := argu[3]; 
		print("den partikulære løsning bliver da");
		if d = sin(argsin) then
			lambda := argsin*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Im(subs(s = lambda/t, H(s))*exp(lambda)));
		elif d = cos(argcos) then
			lambda := argcos*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Re(subs(s = lambda/t, H(s))*exp(lambda)));
		elif d = exp(argexp) then
			lambda := argexp;
#			print("den partikulære løsning bliver da"); 
			y := subs(s = lambda/t, H(s))*exp(lambda); 
		elif d = cos(argcos)*exp(argexp) then
			lambda := argexp + argcos*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Re(subs(s = lambda/t, H(s))*exp(lambda)));
		elif d = sin(argsin)*exp(argexp) then 
			lambda := argexp + argsin*I;
#			print("den partikulære løsning bliver da");
			y := evalc(Im(subs(s = lambda/t, H(s))*exp(lambda))); 
		end if;
	end proc


end module: # Math 2
;
with(math2);
today := Date(timezone = Calendar:-HostTimeZone()):
cat(DayOfMonth(today), ". ", ["januar", "februar", "marts", "april", "maj", "juni", "juli", "august", "september", "oktober", "november", "december"][Month(today)], " ", Year(today));
