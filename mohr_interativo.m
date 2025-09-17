function mohr_interativo()
% Estados
  S.sx    = 80;     % MPa
  S.sy    = -40;    % MPa
  S.txy   = 25;     % MPa  (horário +)
  S.gamma = 0;      % graus
  S.alpha = 0;      % graus
  S.kArrow= 0.50;   % escala das setas
  S.eps   = 1e-12;
  S.dragMode = 'none';  S.theta0 = 0;  S.alpha0 = 0;

  %cores e tamanho
  S.col.outline = [0 0 0];
  S.col.sx  = [0.85 0.00 0.00];
  S.col.sy  = [1.00 0.50 0.00];
  S.col.txy = [0.00 0.55 0.00];
  S.col.pr1 = [0.60 0.00 0.60];
  S.col.pr2 = [0.00 0.45 0.85];
  S.col.tmx = [0.00 0.65 0.80];
  S.ms      = 12;          
  S.fs      = 13;         
  S.fsLab   = 15;        
  S.offLab  = 1.6;         

  %interface
  S.fig = figure('Name','Mohr Interativo','NumberTitle','off', ...
                 'Color','w','Units','normalized','Position',[0.06 0.08 0.88 0.82]);

  % eixos
  S.axElem = subplot(1,2,1,'Parent',S.fig);
  setupAxes(S.axElem,S.fs);
  title(S.axElem,'Elemento (arraste para girar α)  —  τ horário +','FontSize',S.fsLab);

  S.axMohr = subplot(1,2,2,'Parent',S.fig);
  setupAxes(S.axMohr,S.fs);
  xlabel(S.axMohr,'σ','FontSize',S.fsLab); ylabel(S.axMohr,'τ  (horário +)','FontSize',S.fsLab);
  title(S.axMohr,'Círculo de Mohr','FontSize',S.fsLab);

  %  coontroles
  mkText(0.02,0.955,'σx (MPa)',S.fs);     S.edSx   = mkEdit(0.08,0.952,S.sx,@onEditState,S.fs);
  mkText(0.16,0.955,'σy (MPa)',S.fs);     S.edSy   = mkEdit(0.22,0.952,S.sy,@onEditState,S.fs);
  mkText(0.30,0.955,'τxy (MPa, horário +)',S.fs); S.edTxy  = mkEdit(0.48,0.952,S.txy,@onEditState,S.fs);
  mkText(0.66,0.955,'γx (graus, horário +)',S.fs); S.edGam  = mkEdit(0.83,0.952,S.gamma,@onEditGamma,S.fs);

  mkText(0.02,0.915,'α (graus, horário +)',S.fs);
  S.edAlp  = mkEdit(0.14,0.912,S.alpha,@onEditAlpha,S.fs);
  S.slAlpha= uicontrol(S.fig,'Style','slider','Min',-180,'Max',180,'Value',S.alpha, ...
              'Units','normalized','Position',[0.23 0.915 0.32 0.03],'Callback',@onSlideAlpha,'BackgroundColor',[0 0 0]);
  S.btZero = uicontrol(S.fig,'Style','pushbutton','String','α = 0°', ...
              'Units','normalized','Position',[0.57 0.912 0.08 0.035],'Callback',@(~,~)setAlpha(0),'FontSize',S.fs);

  S.txtRead = uicontrol(S.fig,'Style','text','String','', ...
    'Units','normalized','Position',[0.68 0.912 0.30 0.035], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontName','Consolas','FontSize',S.fs);

  % arrastar com mouse
  set(S.fig,'WindowButtonDownFcn',@onMouseDown);
  set(S.fig,'WindowButtonMotionFcn',@onMouseMove);
  set(S.fig,'WindowButtonUpFcn',@onMouseUp);

  guidata(S.fig,S);
  updateAll(S.fig);
end

%utils interface
function setupAxes(ax,fs)
  set(ax,'Color','w'); axis(ax,'equal'); grid(ax,'on'); hold(ax,'on');
  set(ax,'GridColor',[0.85 0.85 0.85],'GridAlpha',0.6,'MinorGridAlpha',0.3);
  set(ax,'XColor','k','YColor','k','FontSize',fs); box(ax,'off');
end
function h = mkText(x,y,str,fs)
  h = uicontrol(gcf,'Style','text','String',str,'Units','normalized', ...
    'Position',[x y 0.15 0.03],'BackgroundColor','w','HorizontalAlignment','left','FontSize',fs);
end
function h = mkEdit(x,y,val,cb,fs)
  h = uicontrol(gcf,'Style','edit','String',num2str(val),'Units','normalized', ...
    'Position',[x y 0.08 0.035],'BackgroundColor','w','Callback',cb,'FontSize',fs);
end

% utils
function onEditState(h,~)
  S = guidata(ancestor(h,'figure'));
  sx = str2double(get(S.edSx,'String')); sy = str2double(get(S.edSy,'String')); tx = str2double(get(S.edTxy,'String'));
  if any(~isfinite([sx,sy,tx])), return; end
  S.sx=sx; S.sy=sy; S.txy=tx; guidata(S.fig,S); updateAll(S.fig);
end
function onEditGamma(h,~)
  S = guidata(ancestor(h,'figure'));
  g = str2double(get(S.edGam,'String')); if ~isfinite(g), return; end
  S.gamma = clamp(g,-180,180); set(S.edGam,'String',sprintf('%.6g',S.gamma));
  guidata(S.fig,S); updateAll(S.fig);
end
function onEditAlpha(h,~)
  S = guidata(ancestor(h,'figure'));
  a = str2double(get(h,'String')); if ~isfinite(a), return; end
  S.alpha = clamp(a,-180,180); set(S.slAlpha,'Value',S.alpha);
  guidata(S.fig,S); updateAll(S.fig);
end
function onSlideAlpha(h,~)
  S = guidata(ancestor(h,'figure'));
  S.alpha = get(h,'Value'); set(S.edAlp,'String',sprintf('%.1f',S.alpha));
  guidata(S.fig,S); updateAll(S.fig);
end
function setAlpha(a)
  S = guidata(gcf);
  S.alpha = clamp(a,-180,180); set(S.edAlp,'String',sprintf('%.1f',S.alpha)); set(S.slAlpha,'Value',S.alpha);
  guidata(S.fig,S); updateAll(S.fig);
end

%arrastar
function onMouseDown(fig,~)
  S = guidata(fig);
  ax = hittest(fig);
  if isempty(ax) || ~ishandle(ax), return; end
  if ancestor(ax,'axes')==S.axElem
      S.dragMode = 'elem';
      cp = get(S.axElem,'CurrentPoint'); p = cp(1,1:2);
      S.theta0 = atan2d(p(2), p(1)); S.alpha0 = S.alpha;
  elseif ancestor(ax,'axes')==S.axMohr
      [C,R] = mohr_center_radius(S.sx,S.sy,S.txy);
      cp = get(S.axMohr,'CurrentPoint'); s = cp(1,1); t = cp(1,2);
      d = hypot(s-C, t-0);
      if abs(d - R) <= 0.20*max(R,1), S.dragMode = 'mohr'; else, S.dragMode = 'none'; end
  else
      S.dragMode = 'none';
  end
  guidata(fig,S);
end
function onMouseMove(fig,~)
  S = guidata(fig);
  switch S.dragMode
    case 'elem'
      cp = get(S.axElem,'CurrentPoint'); p = cp(1,1:2);
      th  = atan2d(p(2), p(1)); dth = wrapDeg(th - S.theta0);
      S.alpha = clamp(S.alpha0 - dth, -180, 180);
      set(S.edAlp,'String',sprintf('%.1f',S.alpha)); set(S.slAlpha,'Value',S.alpha);
      guidata(fig,S); updateAll(fig);
    case 'mohr'
      [C,~] = mohr_center_radius(S.sx,S.sy,S.txy);
      cp = get(S.axMohr,'CurrentPoint'); s = cp(1,1); t = cp(1,2);
      v0 = [S.sx - C, S.txy]; v = [s - C, t];
      if any(~isfinite(v)), return; end
      phi = atan2( v0(1)*v(2) - v0(2)*v(1),  v0(1)*v(1) + v0(2)*v(2) );
      S.alpha = clamp(-0.5*rad2deg(phi),-180,180);
      set(S.edAlp,'String',sprintf('%.1f',S.alpha)); set(S.slAlpha,'Value',S.alpha);
      guidata(fig,S); updateAll(fig);
  end
end
function onMouseUp(fig,~)
  S = guidata(fig); S.dragMode = 'none'; guidata(fig,S);
end

% update tudo
function updateAll(fig)
  S = guidata(fig);
  if any(~isfinite([S.sx,S.sy,S.txy,S.alpha,S.gamma])), return; end

 
  [C,R]       = mohr_center_radius(S.sx,S.sy,S.txy);
  [s1,s2,ap]  = principal_angles(S.sx,S.sy,S.txy);

  [sxp,syp,txyp] = transform_full(S.sx,S.sy,S.txy,S.alpha);
  A = [sxp,  txyp];
  B = [syp, -txyp];

  theta_x = -S.gamma + 90;
  [xP,yP] = pole_point(S.sx,S.sy,S.txy,theta_x);

  drawElement(S.axElem,S,sxp,syp,txyp);
  drawMohr(S.axMohr,S,C,R,s1,s2,[xP yP],ap,A,B);

  set(S.txtRead,'String',sprintf('C=%g  R=%g  σ1=%g  σ2=%g  |τmax|=%g  αp=%g  α=%g', ...
      C,R,s1,s2,abs(R),ap,S.alpha));
end

%contas
function [C,R] = mohr_center_radius(sx,sy,txy)
  C = 0.5*(sx+sy);
  R = hypot(0.5*(sx-sy), txy);
end
function [s1,s2,ap_deg] = principal_angles(sx,sy,txy)
  [C,R] = mohr_center_radius(sx,sy,txy);
  s1 = C + R; s2 = C - R;
  if abs(sx-sy) < 1e-15 && abs(txy) < 1e-15, ap_deg = 0;
  else, ap_deg = 0.5*rad2deg(atan2(2*txy, sx - sy)); % horário +
  end
end
function [sxp,syp,txyp] = transform_full(sx,sy,txy,alpha_deg)
  a = deg2rad(alpha_deg);
  C = 0.5*(sx+sy); D = 0.5*(sx-sy);
  c2 = cos(2*a); s2 = sin(2*a);
  sxp  = C + D*c2 + txy*s2;
  syp  = C - D*c2 - txy*s2;
  txyp = -D*s2 + txy*c2;
end
function [xP,yP] = pole_point(sx,sy,txy,theta_x_deg)
  [C,R] = mohr_center_radius(sx,sy,txy);
  if R <= 1e-9, xP = C+R; yP = 0; return; end
  cBeta = (sy - C)/R; sBeta =  txy / R;
  th = deg2rad(theta_x_deg); c2t = cos(2*th); s2t = sin(2*th);
  xP = C + R*( cBeta*c2t - sBeta*s2t );
  yP =     R*( cBeta*s2t + sBeta*c2t );
end

%design
function drawElement(ax,S,sxp,syp,txyp)
  cla(ax); setupAxes(ax,S.fs);
  title(ax,'Elemento (arraste para girar α)  —  τ horário +','FontSize',S.fsLab);

  % quadrado
  L = 1.0; p = 0.5*L; lim = 1.22;
  P = [ -p -p;  p -p;  p  p; -p  p; -p -p ]';
  phi = deg2rad(-S.alpha); Rz = [cos(phi) -sin(phi); sin(phi) cos(phi)];
  Q = Rz*P; plot(ax,Q(1,:),Q(2,:),'-','Color',S.col.outline,'LineWidth',2);


  magRef = max([abs(sxp),abs(syp),abs(txyp),1]);
  scale  = S.kArrow * L / magRef;

  mid_x = Rz*[ p; 0];   n_x = Rz*[ 1; 0];   t_x = Rz*[ 0;-1];
  mid_y = Rz*[ 0; p];   n_y = Rz*[ 0; 1];   t_y = Rz*[ 1; 0];

  %vetores
  vNx = scale*sxp * n_x;  vNy = scale*syp * n_y;
  vTx = scale*txyp* t_x;  vTy = scale*txyp* t_y;

  quiver(ax, mid_x(1),mid_x(2), vNx(1),vNx(2), 0, 'Color',S.col.sx,'LineWidth',2.2,'MaxHeadSize',0.8);
  quiver(ax, mid_y(1),mid_y(2), vNy(1),vNy(2), 0, 'Color',S.col.sy,'LineWidth',2.2,'MaxHeadSize',0.8);
  quiver(ax, mid_x(1),mid_x(2), vTx(1),vTx(2), 0, 'Color',S.col.txy,'LineWidth',2.2,'MaxHeadSize',0.8);
  quiver(ax, mid_y(1),mid_y(2), vTy(1),vTy(2), 0, 'Color',S.col.txy,'LineWidth',2.2,'MaxHeadSize',0.8);

 
  off = S.offLab;
  text(ax, mid_x(1)+off*vNx(1), mid_x(2)+off*vNx(2), sprintf('σx′=%.3g',sxp), ...
       'FontName','Consolas','Color',S.col.sx,'HorizontalAlignment','center','FontSize',S.fsLab);
  text(ax, mid_y(1)+off*vNy(1), mid_y(2)+off*vNy(2), sprintf('σy′=%.3g',syp), ...
       'FontName','Consolas','Color',S.col.sy,'HorizontalAlignment','center','FontSize',S.fsLab);
  text(ax, mid_x(1)+off*vTx(1), mid_x(2)+off*vTx(2), sprintf('τx′y′=%.3g',txyp), ...
       'FontName','Consolas','Color',S.col.txy,'HorizontalAlignment','center','FontSize',S.fsLab);
  text(ax, mid_y(1)+off*vTy(1), mid_y(2)+off*vTy(2), sprintf('τx′y′=%.3g',txyp), ...
       'FontName','Consolas','Color',S.col.txy,'HorizontalAlignment','center','FontSize',S.fsLab);

  axis(ax,lim*[-1 1 -1 1]);


  hL = [ ...
    plot(ax,nan,nan,'-','Color',S.col.outline,'LineWidth',2), ...
    plot(ax,nan,nan,'-','Color',S.col.sx,'LineWidth',2.2), ...
    plot(ax,nan,nan,'-','Color',S.col.sy,'LineWidth',2.2), ...
    plot(ax,nan,nan,'-','Color',S.col.txy,'LineWidth',2.2) ...
  ];
  legend(ax,hL,{'Borda','σx′','σy′','τx′y′'},'Location','southoutside','Orientation','vertical','Box','off','FontSize',S.fs);
end

%respostas
function drawMohr(ax,S,C,R,s1,s2,P,ap,A,B)
  cla(ax); setupAxes(ax,S.fs);
  xlabel(ax,'σ','FontSize',S.fsLab); ylabel(ax,'τ  (horário +)','FontSize',S.fsLab);
  title(ax,'Círculo de Mohr','FontSize',S.fsLab);

  sigmin = min([s2, A(1), B(1), C-R]); sigmax = max([s1, A(1), B(1), C+R]);
  padX   = 0.15*max(1, sigmax - sigmin);
  padY   = 0.15*max(1, 2*R);
  plot(ax,[sigmin-padX, sigmax+padX],[0,0],'k-','LineWidth',2);
  plot(ax,[0,0],[-R-padY, R+padY],'k-','LineWidth',2);

  th = linspace(-pi,pi,400);
  x = C + R*cos(th); y = R*sin(th);
  plot(ax,x,y,'-','Color',S.col.outline,'LineWidth',2);

  plot(ax,s1,0,'bo','MarkerFaceColor','b','MarkerSize',S.ms); text(ax,s1,0,'  σ1','Color','b','FontSize',S.fsLab);
  plot(ax,s2,0,'bo','MarkerFaceColor','b','MarkerSize',S.ms); text(ax,s2,0,'  σ2','Color','b','FontSize',S.fsLab);
  plot(ax,C,R,'kd','MarkerFaceColor','c','MarkerSize',S.ms);   text(ax,C,R,'  τmax','Color','c','FontSize',S.fsLab);
  plot(ax,C,-R,'kd','MarkerFaceColor','c','MarkerSize',S.ms);  text(ax,C,-R,'  -τmax','Color','c','FontSize',S.fsLab);

  plot(ax,C,0,'ko','MarkerFaceColor','k','MarkerSize',S.ms-2); text(ax,C,0,'  C','Color','k','FontSize',S.fsLab);
  plot(ax,P(1),P(2),'kp','MarkerFaceColor','y','MarkerSize',S.ms+2); text(ax,P(1),P(2),'  P','Color','k','FontSize',S.fsLab);

  plot(ax,A(1),A(2),'ko','MarkerFaceColor','k','MarkerSize',S.ms); text(ax,A(1),A(2),'  A','Color','k','FontSize',S.fsLab);
  plot(ax,B(1),B(2),'ko','MarkerFaceColor','k','MarkerSize',S.ms); text(ax,B(1),B(2),'  B','Color','k','FontSize',S.fsLab);

  %linhas
  plot(ax,[A(1) B(1)],[A(2) B(2)],'-','Color',[0.25 0.25 0.25],'LineWidth',1.8);   % AB
  plot(ax,[C A(1)],[0 A(2)],':','Color',[0.4 0.4 0.4],'LineWidth',1.4);            % C–A
  plot(ax,[C B(1)],[0 B(2)],':','Color',[0.4 0.4 0.4],'LineWidth',1.4);            % C–B
  plot(ax,[P(1) A(1)],[P(2) A(2)],'--','Color',[0.2 0.2 0.2],'LineWidth',1.6);     % P–A
  plot(ax,[P(1) B(1)],[P(2) B(2)],'--','Color',[0.2 0.2 0.2],'LineWidth',1.6);     % P–B
  plot(ax,[P(1) s1],[P(2) 0],'-','Color',S.col.pr1,'LineWidth',2);                  % Psigma1
  plot(ax,[P(1) s2],[P(2) 0],'-','Color',S.col.pr2,'LineWidth',2);                  % PPsigma12
  plot(ax,[P(1) C],[P(2) R], '-','Color',S.col.tmx,'LineWidth',2);                  % P+taumax
  plot(ax,[P(1) C],[P(2) -R],'-','Color',S.col.tmx,'LineWidth',2);                  % P-taymax

  axis(ax,[sigmin-padX, sigmax+padX, -R-padY, R+padY]);

  f = @(x) sprintf('%.3g',x);
  ap2 = wrapDeg(ap+90);
  hAns = [ ...
    plot(ax,nan,nan,'w.'), ...
    plot(ax,nan,nan,'w.'), ...
    plot(ax,nan,nan,'w.'), ...
    plot(ax,nan,nan,'w.') ...
  ];
  legAns = { ...
    ['Planos principais:  αp = ',f(ap),'°  e  αp+90° = ',f(ap2),'°'], ...
    ['Tensões principais:  σ1 = ',f(s1),' MPa,  σ2 = ',f(s2),' MPa'], ...
    ['Cisalhamento máx.:  |τmax| = ',f(abs(R)),' MPa'], ...
    ['Normal correspondente:  σ = C = ',f(C),' MPa'] ...
  };
  lg = legend(ax,hAns,legAns,'Location','southoutside','Orientation','vertical','Box','on','FontSize',S.fsLab);
  lg.ItemTokenSize = [1 1];
end

function v = clamp(v,a,b), v = max(a,min(b,v)); end
function d = wrapDeg(d), d = mod(d+180,360)-180; end
function r = deg2rad(d), r = d*pi/180; end
function d = rad2deg(r), d = r*180/pi; end
