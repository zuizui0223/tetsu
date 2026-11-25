import unittest
from scan_world_railstatus_bird_avonet_guilds import group_pair,MIN_GUILD_HALF
from scan_world_railstatus_observer_split_beta import bucket,add,observer_side
from analyze_world_railstatus_bird_avonet_connectivity import guild_pair
class TestGuild(unittest.TestCase):
    def test_disjoint_guilds(self):
        a=bucket();b=bucket()
        for side in (0,1):
            ids=[i for i in range(100) if observer_side(i)==side][:3]
            for i in range(15):
                add(a,i%7,ids[i%3]);add(b,i%7,ids[i%3])
        z=group_pair(a,b);self.assertEqual(z["split_similarity"]["01"]["jaccard"],1)
        self.assertIsNone(group_pair(a,bucket()))
        r=guild_pair({},{"guilds":{"dense":z,"open":z}},("dense","open"))
        self.assertEqual(r["difference_jaccard"],0.)
if __name__=="__main__":unittest.main()
