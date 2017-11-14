        # why are we doing this?
        if commit:
            commit_consistency_check = commit_id[0:len(commit)]==commit
            if not commit_consistency_check:
                raise ValueError('commit id specified and commit id actually used are not the same:' +commit+commit_id[0:len(commit)])

        self.commit_id = commit_id

        try:
            git_dir = GFDL_BASE+'/.git'
            commit_id_base = sh.git("--git-dir="+git_dir, "log", "--pretty=format:'%H'", "-n 1")
            commit_id_base = str(commit_id_base).split("'")[1]
            git_diff_output = sh.git("--no-pager", "--git-dir="+git_dir, "--work-tree="+self.srcdir, "diff", "--no-color", self.srcdir)
            git_diff_output = str(git_diff_output).split("\n")
        except:
            commit_id_base = None
            git_diff_output  = None

        self.git_diff_output   = git_diff_output

        if commit_id==commit_id_base:
            self.commit_id_base = self.commit_id
        else:
            self.commit_id_base = commit_id_base

        if not (repo and commit):
            try:
                git_dir = self.srcdir+'/.git'
                git_status_output = sh.git("--git-dir="+git_dir, "--work-tree="+self.srcdir, "status", "-b", "--porcelain")
                git_status_output = str(git_status_output)
                git_status_output = git_status_output.split("\n")
                git_status_final = ["Running from GFDL_BASE repo, so adding git status output.\n", git_status_output[0]]
                for suffix_to_search_for in ['.f90', '.inc']:
                    git_status_relevant = [str_entry for str_entry in git_status_output if suffix_to_search_for in str_entry.lower()]
                    git_status_final.extend(git_status_relevant)
            except:
                git_status_final = None
        else:
            git_status_final = ['run using specific commit, as specified above, so no git status output for f90 or inc files']

        self.git_status_output = git_status_final
