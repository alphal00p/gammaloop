#import "shared.typ": *
#import "components.typ": *

#let months = (
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
)

#let date-label(date) = {
  let parts = date.split("-")
  parts.at(2) + " " + months.at(int(parts.at(1)) - 1)
}

#let talk-entry(talk, year) = {
  let person = portal.people.find(person => person.id == talk.speaker)
  talk-card(talk.id, talk.date, date-label(talk.date), year)[
      #kicker[#link("../people/#" + person.id)[#person.name]]

      === #talk.title

      #strong(talk.event) \
      #talk.location

      #navigation("Resources for " + talk.title)[
        #link(talk.event_url)[Event record #arrow()]
        #if "slides_url" in talk {
          link(talk.slides_url)[Slides #arrow()]
        }
        #if "recording_url" in talk {
          link(talk.recording_url)[Recording #arrow()]
        }
      ]
  ]
}

#let sorted-talks = talk-catalog.talk.sorted(
  key: talk => talk.date,
  by: (left, right) => left >= right,
)
#let years = sorted-talks.map(talk => talk.date.slice(0, 4)).dedup()

#let talks = [
  #talks-hero[
    #kicker[Seminars & conferences]

    = Talks

    Selected presentations on Local Unitarity, numerical perturbation theory, GammaLoop, and the scientific software surrounding them.

    #talks-provenance[#str(sorted-talks.len()) talks · linked to public conference records, slides, and recordings where available.]
  ]
  #talk-timeline[
    #for year in years {
      talk-year(year)[
        #for talk in sorted-talks.filter(talk => talk.date.starts-with(year)) {
          talk-entry(talk, year)
        }
      ]
    }
  ]
]

#let talks-page = portal-document(
  "Talks",
  "Talks by αLoop collaborators on Local Unitarity, numerical methods, GammaLoop, and related research software.",
  talks,
  active: "talks",
)
