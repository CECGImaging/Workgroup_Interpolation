function [fwd_bsp,fwd_mat,opt] = Handle_BEM_FWD(torso,heart,EGM,settings)
% Handles the SCIRun boundary element formulation iwth linear mapping over elements using a trisurf mesh.
if isfield(torso,'pts')
    torso.node = torso.pts;
    torso = rmfield(torso,'pts');
end
if isfield(heart,'pts')
    heart.node = heart.pts;
    heart = rmfield(heart,'pts');
end
if isfield(torso,'fac')
    torso.face = torso.fac;
    torso = rmfield(torso,'fac');
end
if isfield(heart,'fac')
    heart.face = heart.fac;
    heart = rmfield(heart,'fac');
end

if ~isfield(settings,'condTorso')
    settings.condTorso = [0,1];
end
if ~isfield(settings,'condHeart')
    settings.condHeart = [1,0];
end


body = struct('surface',[]);
body.surface{1} = struct('node',torso.node,'face',torso.face,'sigma', settings.condTorso);
body.surface{2} = struct('node',heart.node,'face',heart.face,'sigma', settings.condHeart);

tmpFwd = real( bemMatrixPP2(body) );
opt.fullForward = tmpFwd;
fwd_mat = tmpFwd(settings.leadLinks{1},settings.leadLinks{2});

fwd_bsp = fwd_mat*EGM;


end

