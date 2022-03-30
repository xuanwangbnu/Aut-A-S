function mask = genmask(Npre,Npost,con,dir,aut)
	%con is the connection probability
	%cond and ko were stuff i removed, sorry
	mask = rand(Npost,Npre)<con;
	
	%if not directed (ie gap junctions), make symmetrical
	if not(dir) && Npre==Npost
		mask = triu(mask);
		mask = mask + mask.';
		mask = mask - diag(diag(mask));
	end
	
	%aut = 1 if autapses not allowed
	if aut && Npre==Npost
		mask = mask - diag(diag(mask));
	end
    mask=double(mask);
% 	mask = mask';
end
