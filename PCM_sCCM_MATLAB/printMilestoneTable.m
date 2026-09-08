function printMilestoneTable(names, results)
%PRINTMILESTONETABLE  Print the t50/t90/t95/t99 console table used by
%runAll.m for a cell array of scenario names and their matching
%ccmSlipModel/radialEnthalpyModel result structs.
fprintf('%-38s %18s %18s %18s %18s\n', 'scenario', 't50', 't90', 't95', 't99');
for k = 1:numel(names)
    m = results{k}.milestones;
    fprintf('%-38s %18s %18s %18s %18s\n', names{k}, ...
        fmtT(m.t50), fmtT(m.t90), fmtT(m.t95), fmtT(m.t99));
end
end
