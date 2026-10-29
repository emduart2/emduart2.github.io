---
name: update-website
description: Add or change content on Eliana's website (news, papers, talks, students, courses, CV, Come to Porto), preview it, and publish it after her approval.
---

# Update Eliana's website

Follow CLAUDE.md throughout: short messages, explain each step, one challenging question, ask before commit and before publishing.

## 1. Start

1. `git switch master && git pull`, to get any changes made on github.com.
2. Create a branch for the change, e.g. `git switch -c update/<short-name>`.
3. Start the preview: `preview_start` with the `website` configuration (Jekyll on port 4000). If it is missing from `.claude/launch.json`, add it: `bundle exec jekyll serve --port 4000 --livereload` run in this folder.

## 2. Make the change

| Request | What to do |
|---|---|
| News item | Add at the top of `_data/news.yml` (date + Markdown text). |
| New paper | Add to `_data/publications.yml` in date order. Look up the venue and DOI on Crossref (`api.crossref.org/works?query.bibliographic=…`) and the arXiv number. Numbering is automatic. |
| Paper published | Change `type` to `peer-reviewed`, add `venue` and `doi`. |
| Coauthor link | Add `Name: URL` to `_data/coauthors.yml`. The name must match `publications.yml` exactly. Open the page to confirm it is the right person. |
| Talk or event | Add a bullet at the top of the right section in `activities/index.md`. |
| Student | Edit the Students section of `teaching/index.md` (programme, thesis title, co-supervisor). |
| New academic year | Eliana opens SIGARRA in the browser pane and signs in herself. Read (never change) *Distribuição de Serviço* for the new year. Update `current_year` and `current` in `_data/porto_courses.yml`, merge last year's courses into `past` (name once, year + semester), and update the teaching news item. |
| Come to Porto, bio, other prose | Draft the text in chat first. Apply only when she says it is set. |
| Role or title | `_config.yml`, then restart the preview. |

## 3. Check

- Look at the changed page in the preview. After CSS changes, reload with the cache bypassed.
- Check that no "@" email address or `mailto:` link was added, and that the new links load.
- For bigger changes, also check at phone width (`resize_window` mobile), then reset to desktop.

## 4. Publish (each step needs her OK)

1. Commit on the branch, ending the message with the Co-Authored-By line.
2. Merge into `master`. If GitHub's `master` changed meanwhile, `git pull` first and resolve any conflict with her.
3. `git push origin master`. GitHub rebuilds in about 1–2 minutes: watch `gh api repos/emduart2/emduart2.github.io/pages/builds/latest` in the background.
4. Open https://emduart2.github.io to confirm, then delete the branch.

To test a production build locally: `PAGES_REPO_NWO=emduart2/emduart2.github.io JEKYLL_ENV=production bundle exec jekyll build -d <scratch dir>`.
