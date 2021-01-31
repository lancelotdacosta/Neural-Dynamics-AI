function hello = simplex2(NatDP, MDP)


%% First lets draw a 2-simplex (three vertices). 
line_width = 2; 
k=2; %2-simplex
simplex_vertices = eye(k+1);
color = linspecer(3);


hello =0;
for trial = 1:size(NatDP,2)
    for i = 1:size(NatDP(trial).xn{1,1,1,1},3)
        for j =1:size(NatDP(trial).xn{1,1,1,1},4)
            NatDP(trial).xn{1,1,1,1}(:,:,i,j);
            MDP(trial).xn{1,1,1,1}(:,:,i,j);
            temp = abs(cat_inf_length(NatDP(trial).xn{1,1,1,1}(:,:,i,j))-...
                    cat_inf_length(MDP(trial).xn{1,1,1,1}(:,:,i,j)));
            if(temp>hello)
                hello = temp;
            end
        end
    end
end
if(hello>2);
    h=1;
    for trial = 1:size(NatDP,2)
        for i = 1:size(NatDP(trial).xn{1,1,1,1},3)
            for j =1:size(NatDP(trial).xn{1,1,1,1},4)
               % NatDP(trial).xn{1,1,1,1}(:,:,i,j);
                %MDP(trial).xn{1,1,1,1}(:,:,i,j);
                
                nat_updates = NatDP(trial).xn{1,1}(:,:,i,j);
                inf_updates = MDP(trial).xn{1,1}(:,:,i,j);

                dif_l= abs(cat_inf_length(NatDP(trial).xn{1,1,1,1}(:,:,i,j))-...
                    cat_inf_length(MDP(trial).xn{1,1,1,1}(:,:,i,j)));
                dif_init = sum(abs(nat_updates(1,:)-inf_updates(1,:)));
                dif_end = sum(abs(nat_updates(end,:)-inf_updates(end,:)));
                if(dif_l>0.3 ...
                && dif_init <0.05 ...
                && dif_end <0.05)
                disp(trial)
                disp(i)
                disp(j)

                %% for plotting
                figure(h), clf,
                simp_vert = [simplex_vertices, simplex_vertices(:,1)];
                plot3(simp_vert(1,:),simp_vert(2,:),simp_vert(3,:), 'Color', color(1,:));
                hold on
                % 
                % direction = [0 0 1];
                % rotate(hplot,direction,180)

                %% Now let s generate t within some range
                %trial = randperm(size(NatDP,2),1);
    %             i = randperm(3,1);
    %             j = randperm(3,1);

                %nat updates in red
                %C= ones(size(nat_updates(:,1)));
                scatter3(nat_updates(:,1), nat_updates(:,2), nat_updates(:,3), 'DisplayName', 'Nat. Grad');
                %scatter3(nat_updates(:,1), nat_updates(:,2), nat_updates(:,3), 'Color', color(2,:));
                hold on
                % normal updates in yellow
                scatter3(inf_updates(:,1), inf_updates(:,2), inf_updates(:,3),'DisplayName', 'AI');

                %set(gca,'visible','off')
                view(160,20)
                hold off
                
                legend
                txt = {'Trial' num2str(trial),'i' num2str(i), 'j' num2str(j)};
                text(1,1,txt);
                drawnow
                
                h=h+1;
                pause(1);
                end
            end
        end
    end
end
